"""
Select the best-localized gravitational-wave events whose sky localisation
falls inside the LSST footprint for a given timeline.

With no --date/--mjd given, the timeline defaults to "now through the end of
the current MJD day" rather than a full upcoming night.

Help for usage : eventsSelector.py -h
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import astropy_healpix as ah
import healpy as hp
import numpy as np
import pandas as pd
import requests
from astropy.time import Time
from ligo.skymap import postprocess
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree

from GWUtils import models_gw
from GWUtils.define import EVENTS_DIRECTORY, SKYMAP_FITS_DIRECTORY
from GWUtils.models_gw import GWEvent
from GWUtils.query_utils import query_cbc

LOGGER = logging.getLogger("eventsSelector")

OBSLOCTAP_URL = "https://usdf-rsp.slac.stanford.edu/obsloctap/schedule"
CHILE_TZ = ZoneInfo("America/Santiago")
LSST_FOV_RADIUS_DEG = np.sqrt(9.6 / np.pi)  # 9.6 deg^2 circular field of view

# O4b ended 2025-01-28; caps --public-only to GWTC events released through that
# run, so later O4c/O5 superevents (not yet vetted/cataloged) are excluded even
# once GWOSC starts publishing them.
O4B_END_GPS = Time("2025-01-28T16:00:00", scale="utc").gps


# Night window 

def default_night_date(now: datetime | None = None) -> date:
    """Chilean calendar date identifying the relevant LSST night.

    Before local noon we are still inside the night that started the previous
    evening; from noon onward the relevant night is the one starting tonight.
    """
    now = now or datetime.now(tz=CHILE_TZ)
    return now.date() if now.hour >= 12 else now.date() - timedelta(days=1)


def night_start_mjd(night_date: date) -> float:
    """MJD used to query the obsloctap schedule for the given night."""
    midnight = datetime(night_date.year, night_date.month, night_date.day, tzinfo=CHILE_TZ)
    return float(np.floor(Time(midnight).mjd) + 0.75)


def resolve_night(
    date_str: str | None, mjd: float | None = None
) -> tuple[str, float, float | None]:
    """Resolve the --date/--mjd argument to a (timeline label, start MJD, window hours) triple.

    --mjd takes the start MJD directly, bypassing calendar-date parsing
    entirely. Prefer it when coordinating across timezones (e.g. an
    instructor in one timezone assigning nights to participants in
    another): a calendar date is ambiguous without also stating which
    timezone's midnight it refers to (here, Chile/Santiago, since that's
    where LSST observes from), whereas an MJD value is unambiguous.

    With neither --date nor --mjd given, there is no specific night to
    target, so the timeline defaults to "now through the end of the current
    MJD day" (MJD days roll over at UTC midnight) instead of a full
    upcoming night -- the returned window hours reflect just what's left of
    today. The third element of the tuple is non-None only in that default
    case; callers should fall back to --window-hours otherwise.
    """
    if mjd is not None:
        night_date = Time(mjd, format="mjd").to_datetime(timezone=CHILE_TZ).date()
        return f"night of {night_date} (MJD {mjd:.2f})", float(mjd), None
    if date_str is not None:
        night_date = date.fromisoformat(date_str)
        if night_date < default_night_date():
            raise ValueError(f"requested night {night_date} is in the past")
        start_mjd = night_start_mjd(night_date)
        return f"night of {night_date} (MJD {start_mjd:.2f})", start_mjd, None

    start_mjd = float(np.asarray(Time.now().mjd))
    end_mjd = float(np.floor(start_mjd) + 1.0)
    window_hours = (end_mjd - start_mjd) * 24.0
    label = f"current MJD {start_mjd:.4f} through end of MJD day (MJD {end_mjd:.0f})"
    return label, start_mjd, window_hours


def fetch_lsst_pointings(
    start_mjd: float, window_hours: float, executed_only: bool = False
) -> pd.DataFrame:
    """Fetch LSST pointings for the given night from obsloctap.

    obsloctap returns the *planned* schedule, which includes pointings that
    may later be aborted (weather, technical issues) -- each row's
    `execution_status` is "Performed" if it was actually taken, "Aborted"
    otherwise (aborted rows have no real obs_id/timestamps, just placeholder
    values). For a night that hasn't happened yet, nothing is "Performed"
    yet, so the full planned schedule is the right thing to use when picking
    an event to watch in advance -- that's the default here.

    Pass executed_only=True to restrict to pointings actually performed,
    e.g. to check after the fact whether a given (past) night's data could
    plausibly contain an alert for a region you were watching.
    """
    response = requests.get(
        OBSLOCTAP_URL, params={"start": start_mjd, "time": window_hours}, timeout=30
    )
    response.raise_for_status()
    pointings = pd.DataFrame(response.json())
    if pointings.empty:
        raise RuntimeError("obsloctap returned no scheduled pointings for this night")

    n_performed = int((pointings["execution_status"] == "Performed").sum())
    LOGGER.info(
        "%d/%d scheduled pointings were actually executed (execution_status == 'Performed')",
        n_performed, len(pointings),
    )
    if executed_only:
        pointings = pointings[pointings["execution_status"] == "Performed"]
        if pointings.empty:
            raise RuntimeError(
                "no pointings were actually executed this night (it may have been "
                "entirely lost to weather/technical issues) -- rerun without "
                "--executed-only to use the planned schedule instead"
            )
    return pointings


# Footprint geometry

def make_pointing_polygon(
    ra: float, dec: float, radius_deg: float = LSST_FOV_RADIUS_DEG, n: int = 32
) -> Polygon:
    """Approximate a single circular LSST pointing as a polygon."""
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ra_offsets = radius_deg / np.cos(np.radians(dec)) * np.cos(angles)
    dec_offsets = radius_deg * np.sin(angles)
    return Polygon(zip(ra + ra_offsets, dec + dec_offsets))


def build_night_footprint(pointings: pd.DataFrame):
    """Merge all individual pointings of a night into a single footprint geometry."""
    polygons = [make_pointing_polygon(row.s_ra, row.s_dec) for _, row in pointings.iterrows()]
    return unary_union(polygons)


# Candidate events

def is_non_terrestrial(event: GWEvent, terrestrial_threshold: float) -> bool:
    """Basic astrophysical-origin filter on the p_astro classification."""
    if event.classification is None:
        return True
    return (event.classification.terrestrial or 0.0) < terrestrial_threshold


def is_public_gwtc_event(event: GWEvent, cache: bool = False, save_results: bool = True) -> bool:
    """True if `event` matches a GWOSC-released GWTC event through O4b.

    GraceDB superevents don't reliably carry their GWTC name (`gw_id` is
    often unset even for cataloged events), so the match is done by GPS time
    against GWOSC's public event list instead, via `event.resolve_is_public()`
    (see GWUtils.models_gw.GWEvent), which caches the result on `is_public`.

    A cached `is_public` (e.g. from a prior run) skips the GWOSC query
    entirely. With `cache=True` and no cached result yet, the event is
    treated as unknown and excluded (with a warning) rather than querying
    GWOSC live, since cache mode is meant to stay offline -- run once
    without --cache to resolve and cache it. Otherwise, a freshly resolved
    result is persisted immediately via `event.save()` unless
    `save_results=False`, so future --cache runs over the same event don't
    need to re-query GWOSC.
    """
    if event.t_0 is None:
        return False
    gps = Time(event.t_0).gps
    if gps > O4B_END_GPS:
        return False
    if event.is_public is None:
        if cache:
            LOGGER.warning(
                "%s has no cached public/private status yet; excluding from "
                "--public-only (re-run once without --cache to resolve and cache it)",
                event.superevent_id,
            )
            return False
        event.resolve_is_public()
        if save_results:
            event.save()
    return bool(event.is_public)


def load_cached_events(events_dir: Path) -> list[GWEvent]:
    """Load GWEvent records previously cached to disk by GWEvent.save() (see GWUtils.models_gw)."""
    events = []
    for path in sorted(Path(events_dir).glob("*.json")):
        try:
            events.append(GWEvent.model_validate(json.loads(path.read_text())))
        except Exception as e:
            LOGGER.warning("Skipping cached event %s: could not load (%s)", path, e)
    return events


def fetch_candidate_events(
    far_threshold: float,
    terrestrial_threshold: float,
    cache: bool = False,
    events_dir: Path = EVENTS_DIRECTORY,
    public_only: bool = False,
    save_results: bool = True,
) -> list[GWEvent]:
    """Get significant CBC superevents with a ready skymap.

    By default queries GraceDB live. If `cache` is set, events are instead loaded
    from the local JSON cache written by GWEvent.save(); their skymaps are
    expected to already be cached too (see GWUtils.models_gw.has_dl_skymap),
    since cache mode does not hit the network to download missing ones.

    If `public_only` is set, candidates are further restricted to events
    already released in the public GWTC catalog through O4b (see
    `is_public_gwtc_event`), excluding non-public/unvetted superevents. That
    check is cached per-event and, unless `save_results` is False, persisted
    immediately so a later --cache run over the same events skips GWOSC
    entirely.
    """
    if cache:
        events = load_cached_events(events_dir)
        events = [ev for ev in events if ev.skymap_ready and ev.far is not None and ev.far < far_threshold]
        cached, missing_skymap = [], 0
        for ev in events:
            if ev.has_dl_skymap():
                cached.append(ev)
            else:
                missing_skymap += 1
        if missing_skymap:
            LOGGER.warning("Skipping %d cached event(s) with no cached skymap", missing_skymap)
        events = cached
    else:
        events = query_cbc(f"far < {far_threshold} label: SKYMAP_READY", enrich=False)
    events = [ev for ev in events if is_non_terrestrial(ev, terrestrial_threshold)]
    if public_only:
        events = [
            ev for ev in events
            if is_public_gwtc_event(ev, cache=cache, save_results=save_results)
        ]
    return events


# Footprint coverage 

def assess_event_coverage(event: GWEvent, footprint, threshold: float) -> dict | None:
    """Compute how much of an event's 90% credible region falls in the footprint."""
    try:
        skymap, _ = event.load_skymap()
    except Exception as e:
        LOGGER.warning("Skipping %s: could not load skymap (%s)", event.superevent_id, e)
        return None

    nside = ah.npix_to_nside(len(skymap))
    credible_levels = postprocess.find_greedy_credible_levels(skymap)
    pixels_90 = np.where(credible_levels <= 0.90)[0]

    theta, phi = hp.pix2ang(nside, pixels_90)
    ra, dec = np.degrees(phi), 90.0 - np.degrees(theta)

    tree = STRtree([Point(r, d) for r, d in zip(ra, dec)])
    in_fp = np.zeros(len(pixels_90), dtype=bool)
    in_fp[tree.query(footprint, predicate="contains")] = True

    pixel_area = ah.nside_to_pixel_area(nside).to_value("deg2")
    prob_covered = float(skymap[pixels_90][in_fp].sum())

    return {
        "event_id": event.gw_id or event.superevent_id,
        "superevent_id": event.superevent_id,
        "far": event.far,
        "classification": event.classification.most_probable() if event.classification else None,
        "prob_covered": prob_covered,
        "area_90_deg2": len(pixels_90) * pixel_area,
        "covered_deg2": in_fp.sum() * pixel_area,
        "in_footprint": prob_covered >= threshold,
    }


def rank_events(
    events: list[GWEvent], footprint, threshold: float, n_events: int, workers: int
) -> pd.DataFrame:
    """Score every candidate event and return the top N ranked by localisation."""
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(assess_event_coverage, ev, footprint, threshold): ev for ev in events}
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                rows.append(result)

    if not rows:
        return pd.DataFrame()

    covered = pd.DataFrame(rows)
    covered = covered[covered["in_footprint"]]
    return covered.sort_values(
        ["area_90_deg2", "covered_deg2"], ascending=[True, False]
    ).head(n_events).reset_index(drop=True)


# Plotting

def plot_night_overview(
    footprint, ranking: pd.DataFrame, events_by_id: dict, output_dir: Path,
    timeline_label: str | None = None,
) -> Path:
    """Save a Mollweide plot of the footprint with the selected events overlaid."""
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    def to_mollweide(ra_deg, dec_deg):
        return np.radians(180.0 - ra_deg), np.radians(dec_deg)

    fig, ax = plt.subplots(figsize=(14, 7), subplot_kw={"projection": "mollweide"})
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#0d1117")

    geoms = footprint.geoms if hasattr(footprint, "geoms") else [footprint]
    for geom in geoms:
        ra, dec = np.array(geom.exterior.coords).T
        ra_r, dec_r = to_mollweide(ra, dec)
        ax.fill(ra_r, dec_r, color="deepskyblue", alpha=0.25, zorder=2)
        ax.plot(ra_r, dec_r, color="deepskyblue", lw=0.6, zorder=3)

    colors = cm.plasma(np.linspace(0.15, 0.95, max(len(ranking), 1)))
    legend_handles = [Patch(facecolor="deepskyblue", alpha=0.5, label="LSST footprint")]

    for color, row in zip(colors, ranking.itertuples()):
        skymap, _ = events_by_id[row.superevent_id].load_skymap()
        nside = ah.npix_to_nside(len(skymap))
        levels = postprocess.find_greedy_credible_levels(skymap)
        pixels = np.where(levels <= 0.90)[0]
        theta, phi = hp.pix2ang(nside, pixels)
        ra_r, dec_r = to_mollweide(np.degrees(phi), 90.0 - np.degrees(theta))
        ax.scatter(ra_r, dec_r, s=0.5, color=color, alpha=0.4, zorder=4, rasterized=True)
        legend_handles.append(
            Line2D([0], [0], marker="*", color=color, linestyle="none", markersize=8, label=row.event_id)
        )

    def ra_hours_formatter(x, pos):
        ra_deg = (180.0 - np.degrees(x)) % 360.0
        return f"{ra_deg / 15.0:.0f}h"

    ax.xaxis.set_major_formatter(plt.FuncFormatter(ra_hours_formatter))

    ax.grid(True, color="white", alpha=0.15, lw=0.5)
    ax.tick_params(colors="white")
    title = f"GW events in LSST footprint (n={len(ranking)})"
    if timeline_label:
        title += f"\n{timeline_label}"
    ax.set_title(title, color="white", fontsize=14, pad=15)
    ax.legend(
        handles=legend_handles, loc="lower left", framealpha=0.2, labelcolor="white",
        fontsize=8, facecolor="#0d1117", edgecolor="white",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    fig_path = output_dir / "gw_footprint_overlap.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig_path


# CLI 

def configure_skymap_directory(path: Path) -> None:
    """Point GWUtils' skymap cache at a custom directory."""
    path.mkdir(parents=True, exist_ok=True)
    models_gw.SKYMAP_FITS_DIRECTORY = path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='''Script to select existing GW events which sky localisation are to be covered by LSST.\n
        Procedure :\n
        1) Parses the LSST footprint scheduled for the indicated date (has to be in future or current night ;
        defaults to current or next night, depending on Chile time.)\n
        2) Fetches all available GW events and selects the best localized events (defaults to 5 events) matching basic criteria (non-terrestrial).
        Events are hierarchized first by descending total localisation area, then by coverage in the night footprint.\n
        3) Returns the list of events identificators (either GWTC name or superevent ID). Also plots a couple figures for visualisation
        ''',
        epilog="""Designed for the need of 2026 Sao Paulo workshop.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--date", type=str, default=None,
        help="ISO date (YYYY-MM-DD) of the observing night to analyze, interpreted as a "
             "Chile/Santiago calendar date (where LSST observes from). If neither --date "
             "nor --mjd is given, the timeline instead defaults to the current MJD "
             "through the end of the current MJD day (no fixed night window). If you "
             "were assigned a night across timezones (e.g. a Sao Paulo workshop), prefer "
             "--mjd instead to avoid date ambiguity.",
    )
    parser.add_argument(
        "--mjd", type=float, default=None,
        help="Start MJD of the observing night to analyze, used as-is (no timezone "
             "involved). Takes precedence over --date if both are given; this is the "
             "unambiguous way to specify a night assigned to you by someone in a "
             "different timezone.",
    )
    parser.add_argument(
        "--window-hours", type=float, default=None,
        help="Duration in hours of the LSST schedule window fetched from the timeline "
             "start. Defaults to 24h when --date or --mjd is given, or to however many "
             "hours remain until the end of the current MJD day when neither is given.",
    )
    parser.add_argument(
        "--executed-only", action="store_true",
        help="Restrict the footprint to pointings LSST actually executed, skipping "
             "planned-but-aborted ones (e.g. lost to weather/technical issues). Only "
             "meaningful for a night that has already happened -- use this to check "
             "after the fact whether a watched region could plausibly have produced "
             "any alerts. A 'N/M scheduled pointings were actually executed' line is "
             "always logged regardless of this flag.",
    )
    parser.add_argument(
        "--n-events", type=int, default=5,
        help="Number of best-localized events to select.",
    )
    parser.add_argument(
        "--far-threshold", type=float, default=1e-7,
        help="GraceDB false-alarm-rate threshold [Hz] used to query candidate superevents.",
    )
    parser.add_argument(
        "--terrestrial-threshold", type=float, default=0.5,
        help="Maximum terrestrial probability allowed for an event to be kept as astrophysical.",
    )
    parser.add_argument(
        "--coverage-threshold", type=float, default=0.5,
        help="Minimum fraction of an event's 90%% credible probability inside the footprint "
             "for it to be flagged as covered by the LSST night.",
    )
    parser.add_argument(
        "--skymap-dir", type=Path, default=SKYMAP_FITS_DIRECTORY,
        help="Directory used to cache downloaded skymap FITS files.",
    )
    parser.add_argument(
        "--events-dir", type=Path, default=EVENTS_DIRECTORY,
        help="Directory used to cache selected GWEvent JSON records.",
    )
    parser.add_argument(
        "--cache", action="store_true",
        help="Do not query GraceDB; select candidate events from the local JSON/skymap "
             "cache instead (see --events-dir / --skymap-dir). The LSST footprint is "
             "always fetched live regardless of this flag.",
    )
    parser.add_argument(
        "--public-only", action="store_true",
        help="Restrict candidates to events already released in the public GWTC "
             "catalog through O4b, excluding non-public/unvetted superevents. This "
             "always needs network access too (queries GWOSC), unless every "
             "candidate's public/private status was already cached by a prior run "
             "(see --no-save-events).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path.cwd() / 'output' / "lsst_gw_selection",
        help="Directory where the figures and the selected-event summary are written.",
    )
    parser.add_argument(
        "--workers", type=int, default=min(os.cpu_count() or 4, 8),
        help="Number of worker threads used to assess skymap coverage in parallel.",
    )
    parser.add_argument(
        "--no-plots", action="store_true",
        help="Skip generating visualisation figures.",
    )
    parser.add_argument(
        "--no-save-events", action="store_true",
        help="Do not cache the selected GWEvent records to disk, and (with "
             "--public-only) do not persist each candidate's freshly resolved "
             "public/private status either.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    return parser


def main(argv: list[str] | None = None) -> list[str]:
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    LOGGER.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    try:
        timeline_label, start_mjd, default_window_hours = resolve_night(args.date, args.mjd)
    except ValueError as e:
        sys.exit(f"error: {e}")
    window_hours = args.window_hours if args.window_hours is not None else (
        default_window_hours if default_window_hours is not None else 24.0
    )

    LOGGER.info("Selecting events for %s", timeline_label)
    configure_skymap_directory(args.skymap_dir)

    try:
        pointings = fetch_lsst_pointings(start_mjd, window_hours, args.executed_only)
    except (requests.RequestException, RuntimeError) as e:
        sys.exit(f"error: could not fetch LSST schedule: {e}")
    LOGGER.info("Fetched %d scheduled pointings", len(pointings))
    footprint = build_night_footprint(pointings)

    candidates = fetch_candidate_events(
        args.far_threshold, args.terrestrial_threshold, cache=args.cache,
        events_dir=args.events_dir, public_only=args.public_only,
        save_results=not args.no_save_events,
    )
    LOGGER.info("Found %d candidate events", len(candidates))
    events_by_id = {ev.superevent_id: ev for ev in candidates}

    ranking = rank_events(candidates, footprint, args.coverage_threshold, args.n_events, args.workers)
    if ranking.empty:
        LOGGER.warning("No candidate event is covered by this night's footprint")
        return []

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(args.output_dir / "selected_events.csv", index=False)

    if not args.no_save_events:
        for sid in ranking["superevent_id"]:
            events_by_id[sid].save()

    if not args.no_plots:
        fig_path = plot_night_overview(footprint, ranking, events_by_id, args.output_dir, timeline_label)
        LOGGER.info("Saved overview figure to %s", fig_path)

    selected_ids = ranking["event_id"].tolist()
    LOGGER.info("Selected events: %s", selected_ids)
    return selected_ids


def cli() -> None:
    """Console-script entry point: discards main()'s return value so a
    successful run doesn't get reported as a failure (setuptools wraps
    console_scripts targets as `sys.exit(target())`, and `sys.exit` treats
    any non-None, non-int return value as an error)."""
    main()


if __name__ == "__main__":
    main()
