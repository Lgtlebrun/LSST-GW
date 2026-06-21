"""
Show a summary, skymap plot and public web links for a single GW event.

Help for usage : eventVisualisation.py -h
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path

from ligo.gracedb.rest import HTTPError

from GWUtils import models_gw
from GWUtils.define import EVENTS_DIRECTORY, SKYMAP_FITS_DIRECTORY
from GWUtils.models_gw import CBCClassification, GWEvent, GWTCEvent, UncertainQuantity
from GWUtils.query_utils import query_cbc

LOGGER = logging.getLogger("eventVisualisation")

GRACEDB_SUPEREVENT_URL = "https://gracedb.ligo.org/superevents/{sid}/view/"
GWOSC_EVENT_URL = "https://gwosc.org/events/{name}/"


# Fetching 

def configure_skymap_directory(path: Path) -> None:
    """Point GWUtils' skymap cache at a custom directory."""
    path.mkdir(parents=True, exist_ok=True)
    models_gw.SKYMAP_FITS_DIRECTORY = path


def load_cached_event(event_id: str, events_dir: Path) -> GWEvent:
    """Look up a single GWEvent previously cached to disk by GWEvent.save(),
    matching by superevent ID or GWTC name (any version)."""
    target = re.sub(r"-v\d+$", "", event_id).upper()
    for path in sorted(Path(events_dir).glob("*.json")):
        try:
            ev = GWEvent.model_validate(json.loads(path.read_text()))
        except Exception as e:
            LOGGER.warning("Skipping cached event %s: could not load (%s)", path, e)
            continue
        candidates = {ev.superevent_id, ev.gw_id, re.sub(r"-v\d+$", "", ev.gw_id or "")}
        if target in {c.upper() for c in candidates if c}:
            return ev
    raise ValueError(f"no cached event matching {event_id!r} found in {events_dir}")


def fetch_event(
    event_id: str, enrich: bool, cache: bool = False, events_dir: Path = EVENTS_DIRECTORY
) -> GWEvent:
    """Resolve a GraceDB superevent ID (e.g. S191117j) or a GWTC name (e.g. GW170817).

    If `cache` is set, the event is loaded from the local JSON cache written by
    GWEvent.save() instead of querying GraceDB/GWOSC; its skymap is expected
    to already be cached too (see GWUtils.models_gw.has_dl_skymap).
    """
    if cache:
        return load_cached_event(event_id, events_dir)

    if event_id.upper().startswith("GW"):
        return GWTCEvent(event_id)

    matches = query_cbc(event_id, enrich=enrich)
    if not matches:
        raise ValueError(f"no CBC superevent found for {event_id!r}")
    if len(matches) > 1:
        raise ValueError(f"query {event_id!r} matched {len(matches)} events, expected one")
    return matches[0]


# Summary 

def format_quantity(uq: UncertainQuantity | None, default_unit: str = "") -> str | None:
    """Render an UncertainQuantity as 'value (lower/upper) unit'."""
    if uq is None:
        return None
    text = f"{uq.value:.3g}"
    if uq.lower is not None and uq.upper is not None:
        text += f" ({uq.lower:+.2g}/{uq.upper:+.2g})"
    unit = uq.unit or default_unit
    return f"{text} {unit}" if unit else text


def format_classification(classification: CBCClassification) -> str:
    """Render the most probable class plus the full probability breakdown."""
    label, prob = classification.most_probable()
    probs = {
        "BBH": classification.bbh,
        "BNS": classification.bns,
        "NSBH": classification.nsbh,
        "Terrestrial": classification.terrestrial,
    }
    detail = ", ".join(f"{k}={v:.1%}" for k, v in probs.items() if v is not None)
    return f"{label} ({prob:.1%}) [{detail}]" if detail else f"{label} ({prob:.1%})"


def summarize_event(event: GWEvent) -> str:
    """Build a human-readable summary of an event's main characteristics."""
    rows: list[tuple[str, str]] = []

    def add(label: str, value) -> None:
        if value:
            rows.append((label, str(value)))

    add("Superevent ID", event.superevent_id)
    add("GWTC name", event.gw_id)
    add("Catalog", event.catalog)
    add("Coalescence time (UTC)", event.t_0.strftime("%Y-%m-%d %H:%M:%S UTC") if event.t_0 else None)
    add("FAR [Hz]", f"{event.far:.2e}" if event.far else None)
    add("Network SNR", event.network_snr)
    add("Detectors", ", ".join(str(d) for d in event.detectors) if event.detectors else None)
    if event.classification is not None:
        add("Classification", format_classification(event.classification))
    add("Mass 1", format_quantity(event.mass_1, "Msun"))
    add("Mass 2", format_quantity(event.mass_2, "Msun"))
    add("Chirp mass", format_quantity(event.chirp_mass, "Msun"))
    add("Total mass", format_quantity(event.total_mass, "Msun"))
    add("Final mass", format_quantity(event.final_mass, "Msun"))
    add("Luminosity distance", format_quantity(event.luminosity_distance, "Mpc"))
    add("Redshift", format_quantity(event.redshift))
    add("Effective spin (chi_eff)", format_quantity(event.chi_eff))

    width = max(len(label) for label, _ in rows)
    return "\n".join(f"{label:<{width}} : {value}" for label, value in rows)


def event_links(event: GWEvent) -> dict[str, str]:
    """Public web links for this event, where applicable."""
    links = {}
    if event.superevent_id:
        links["GraceDB"] = GRACEDB_SUPEREVENT_URL.format(sid=event.superevent_id)
    if event.gw_id:
        name = re.sub(r"-v\d+$", "", event.gw_id)
        links["GWTC"] = GWOSC_EVENT_URL.format(name=name)
    return links


# Plotting

def plot_skymap(event: GWEvent, output_dir: Path, roi: str, n_vertices: int) -> Path:
    """Render the event's skymap via GWEvent.plot_event, with an optional ROI overlay."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_path = output_dir / f"{event.identifier}_skymap.png"
    event.plot_event(
        figPath=fig_path,
        # circle_roi=roi in ("circle", "both"),
        # rect_roi=roi in ("rect", "both"),
        n_vertices=n_vertices if roi in ("moc", "both") else None,
    )
    return fig_path


def save_roi_moc(event: GWEvent, output_dir: Path, n_vertices: int, percentile: float = 90) -> Path:
    """Save the event's credible region (see GWEvent.get_roi) as a MOC FITS file."""
    moc = event.get_roi(percentile=percentile, n_vertices=n_vertices, format="moc")
    output_dir.mkdir(parents=True, exist_ok=True)
    moc_path = output_dir / f"{event.identifier}_roi.fits"
    moc.save(moc_path, format="fits", overwrite=True)
    return moc_path


# CLI

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='''Script to fetch statistics and visualise localisation maps of an event, indicated by either GraceDB superevent ID (ex : S191117j) or GWTC ID (ex : GW170817) :\n
       ''',
        epilog="""Designed for the need of 2026 Sao Paulo workshop.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "event_id", type=str,
        help="GraceDB superevent ID (e.g. S191117j) or GWTC catalog name (e.g. GW170817).",
    )
    parser.add_argument(
        "--roi", choices=["circle", "rect", "both", "moc", "none"], default="moc",
        help="90%% credible-region overlay drawn on the skymap.",
    )
    parser.add_argument(
        "--n-vertices", type=int, default=50,
        help="Boundary vertex budget for the MOC region of interest (see GWEvent.get_roi), "
             "used for both the optional 'moc' skymap overlay and the saved MOC file.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path.cwd() / 'output' /  "event_visualisation",
        help="Directory where the skymap figure and MOC region file are written.",
    )
    parser.add_argument(
        "--skymap-dir", type=Path, default=SKYMAP_FITS_DIRECTORY,
        help="Directory used to cache downloaded skymap FITS files.",
    )
    parser.add_argument(
        "--events-dir", type=Path, default=EVENTS_DIRECTORY,
        help="Directory used to cache GWEvent JSON records (see --cache).",
    )
    parser.add_argument(
        "--cache", action="store_true",
        help="Do not query GraceDB/GWOSC; look up the event from the local JSON/skymap "
             "cache instead (see --events-dir / --skymap-dir).",
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip generating the skymap figure.",
    )
    parser.add_argument(
        "--no-enrich", action="store_true",
        help="Skip the (slower) GWOSC catalog lookup when given a GraceDB superevent ID.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    return parser


def main(argv: list[str] | None = None) -> GWEvent:
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    LOGGER.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    configure_skymap_directory(args.skymap_dir)

    try:
        event = fetch_event(
            args.event_id, enrich=not args.no_enrich, cache=args.cache, events_dir=args.events_dir
        )
    except (ValueError, HTTPError) as e:
        sys.exit(f"error: {e}")

    print(summarize_event(event))

    links = event_links(event)
    if links:
        print()
        for name, url in links.items():
            print(f"{name}: {url}")

    if not args.no_plot:
        fig_path = plot_skymap(event, args.output_dir, args.roi, args.n_vertices)
        LOGGER.info("Saved skymap figure to %s", fig_path)

        moc_path = save_roi_moc(event, args.output_dir, args.n_vertices)
        LOGGER.info("Saved ROI MOC to %s", moc_path)

    return event


def cli() -> None:
    """Console-script entry point: discards main()'s return value so a
    successful run doesn't get reported as a failure (setuptools wraps
    console_scripts targets as `sys.exit(target())`, and `sys.exit` treats
    any non-None, non-int return value as an error)."""
    main()


if __name__ == "__main__":
    main()
