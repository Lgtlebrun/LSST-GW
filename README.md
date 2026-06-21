# LSST-GW

Tools to cross-match LIGO/Virgo/KAGRA gravitational-wave events with the LSST
observing schedule: select the best-localized GW events visible to LSST on a
given night, and visualise/summarize any individual event.

Built for the 2026 São Paulo LSST workshop.

This package depends on [GWUtils](https://github.com/Lgtlebrun/GWUtils), a
separate helper library for querying GraceDB/GWOSC and handling GW event data.
You don't need to install it yourself — `pip install .` (see below) fetches
and installs it automatically straight from its GitHub repository.

## 1. Prerequisites

- Python 3.10+ (tested with 3.12)
- `git` (needed by pip to fetch the GWUtils dependency)
- A LIGO.ORG account is only required for non-public/embargoed GraceDB events.
  Public events (the default for most use cases) work without any
  credentials.

## 2. Get the code

```bash
git clone <this-repository-url> LSST-GW
cd LSST-GW
```

(If you already have this folder, just `cd` into it.)

## 3. Create and activate a virtual environment

A virtual environment ("venv") keeps this project's Python packages separate
from the rest of your system.

```bash
python3 -m venv .venv
source .venv/bin/activate      # on Windows: .venv\Scripts\activate
```

You'll know it worked because your terminal prompt now starts with `(.venv)`.
Run `source .venv/bin/activate` again any time you open a new terminal and
want to use this project.

## 4. Install the package

```bash
pip install --upgrade pip
pip install .
```

This installs LSST-GW itself, GWUtils (pulled straight from its GitHub repo),
and every other dependency (astropy, healpy, matplotlib, ligo.skymap, etc.).
It can take a few minutes the first time, since some of these packages are
large. If you plan on editing the code, install with `pip install -e .`
instead so your changes take effect without reinstalling.

To check everything installed correctly:

```bash
python -c "import LSST_GW, GWUtils; print('OK')"
```

## 5. Run the scripts

Installing the package also installs two commands directly on your `PATH`:

### `lsst-gw-select-events`

Selects the best-localized GW events whose sky localisation overlaps
tomorrow night's LSST footprint.

```bash
lsst-gw-select-events
```

Useful options:
- `--date YYYY-MM-DD` — analyze a specific night, interpreted as a
  Chile/Santiago calendar date (where LSST observes from). Defaults to the
  next LSST night.
- `--mjd MJD` — analyze a specific night by its start MJD instead, used
  as-is with no timezone conversion. Prefer this over `--date` if you're
  coordinating across timezones (e.g. assigned a night by someone elsewhere)
  to avoid date-boundary ambiguity.
- `--n-events N` — how many events to select (default 5).
- `--public-only` — restrict to events already publicly released in the GWTC
  catalog (through O4b), excluding unvetted/non-public superevents.
- `--executed-only` — restrict the footprint to pointings LSST actually
  executed, excluding planned-but-aborted ones (weather, technical issues).
  Only meaningful for a night that already happened; use it to check after
  the fact whether a watched region could plausibly have produced any
  alerts. The log always reports "N/M scheduled pointings were actually
  executed" regardless of this flag, even without setting it.
- `--cache` — skip GraceDB/GWOSC and use the local JSON/skymap cache instead
  for GW events. This script still always fetches the LSST footprint live
  from the obsloctap schedule, so an internet connection is required even
  with `--cache`.
- `-v` — verbose/debug logging.

Run `lsst-gw-select-events --help` for the full list. Results (a CSV summary
and an overview figure) are written to `./output/lsst_gw_selection/` by
default.

### `lsst-gw-visualise-event`

Prints a summary and plots the skymap for a single event, given either a
GraceDB superevent ID (e.g. `S190503bf`) or a GWTC catalog name (e.g.
`GW170817`).

```bash
lsst-gw-visualise-event GW170817
```

Useful options:
- `--cache` — look up the event from the local cache instead of querying GraceDB/GWOSC.
- `--no-plot` — skip generating the skymap figure.
- `-v` — verbose/debug logging.

Run `lsst-gw-visualise-event --help` for the full list. Output (skymap figure
and region-of-interest file) is written to `./output/event_visualisation/` by
default.

## Troubleshooting

- **`pip install .` fails while building `healpy`/`lalsuite`/`ligo.skymap`**:
  these packages ship pre-built wheels for common platforms (Linux, macOS on
  Python 3.10-3.12); if pip falls back to building from source, make sure you
  have a C/Fortran compiler available, or try a supported Python version.
- **`401 Unauthorized` errors when fetching an event**: that event isn't
  public yet and requires LIGO.ORG credentials (a valid scitoken or X.509
  certificate) configured on your machine. Public, already-released events
  (most of GWTC) don't need this.
- **A live run seems slow**: both scripts query GraceDB live by default,
  which can be slow if many candidate superevents match (each one needs a
  follow-up request for its classification). Narrow the search with
  `--far-threshold` (for `lsst-gw-select-events`) or use `--cache` once
  you've cached the events you need.
- **`lsst-gw-select-events` hangs or times out fetching the schedule**: this
  command always fetches the LSST footprint live from the obsloctap service,
  even with `--cache` (see above), so it needs a working internet
  connection. On a network that silently drops outbound connections (e.g. a
  restrictive firewall) the request can take noticeably longer than its
  30s timeout before failing — if it never returns, check that you can
  reach `usdf-rsp.slac.stanford.edu`.
- **No alerts came in overnight for a region you watched**: `lsst-gw-select-events`
  builds the footprint from the *planned* LSST schedule, not from what was
  actually observed — some planned pointings get aborted (weather,
  technical issues), and on a bad night that can be the entire schedule.
  Re-run the selector for that night with `--executed-only` to check
  whether there was any real coverage at all, independent of which region
  you were watching.
