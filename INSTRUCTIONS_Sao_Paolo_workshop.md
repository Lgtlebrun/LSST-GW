# São Paulo LSST-GW Workshop — Instructions

This workshop walks through the full chain: pick a real gravitational-wave
(GW) event, build a region-of-interest (watchmap) for it, configure a Rubin
broker to watch that region, and the next day check what LSST actually saw
there.

## Day 1

### 1. Install the package

Follow [README.md](README.md) (clone, create a venv, `pip install .`). At
the end, check the install with:

```bash
python -c "import LSST_GW, GWUtils; print('OK')"
```

If `pip install .` fails on `healpy`/`lalsuite`/`ligo.skymap`, or you hit a
`401 Unauthorized` fetching an event, see the Troubleshooting section at the
bottom of the README.

### 2. Select an event

> The professor will give live instructions on group composition and which
> event/night is assigned to your group. **Double-check you select the
> correct night** — `lsst-gw-select-events` defaults to the *current or
> next* LSST night (Chile time), which may not be the one assigned to you.

Run the selector for your assigned night:

```bash
lsst-gw-select-events --date YYYY-MM-DD
```

This fetches the LSST schedule for that night, ranks GW superevents whose
sky localisation overlaps the night's footprint, and writes:
- `output/lsst_gw_selection/selected_events.csv` — ranked candidate events
- `output/lsst_gw_selection/gw_footprint_overlap.png` — overview plot

Useful flags: `--n-events N`, `--public-only` (restrict to GWTC-cataloged
events since you are most likely not a LIGO member and are thus not allowed to fetch non-public events, which would raise a 401 error), `--offline` (use the local cache instead of querying
GraceDB/GWOSC live, which is necessaron first run, but slow afterwards). Run `lsst-gw-select-events --help` for the full list.

Pick one event ID (superevent ID, e.g. `S190503bf`, or GWTC name, e.g.
`GW170817`) from the CSV to carry forward to the next step.

### 3. Visualise the event and generate a MOC

```bash
lsst-gw-visualise-event <event_id>
```

This prints a summary (masses, distance, FAR, classification...), saves the
skymap plot, and saves a Multi-Order Coverage (MOC) FITS file of the 90%
credible region — this MOC file is what you'll upload to Lasair in step 4.
Output goes to `output/event_visualisation/` by default.

Useful flags:
- `--roi {circle,rect,both,moc,none}` — overlay style on the skymap plot ; the one you need is `moc`
- `--n-vertices N` — boundary complexity of the MOC region (default 50)
- `--offline` — use the local cache instead of querying GraceDB/GWOSC

For a more hands-on, step-by-step walkthrough of the same operations (load
an event, plot it, extract the MOC as a `dict`/`mocpy.MOC`, save it), open
[GWtutorial.ipynb](GWtutorial.ipynb) and play around with it. In particular:

```python
from GWUtils.models_gw import GWEvent

event = GWEvent("<event_id>", offline=False)
event.plot_event(n_vertices=50)

# Inspect the event's full specs
event.to_dict()        # flat dict, ready for inspection or to_dataframe()
event.to_dataframe()    # one-row pandas DataFrame
```

`event.get_roi(n_vertices=50, format="moc")` returns the `mocpy.MOC` object
that you save with `moc.save(path, format="fits", overwrite=True)`; this is
the same file `lsst-gw-visualise-event` writes for you automatically.

### 4. Configure Lasair for the night

Follow [Lasair tutorial.pdf](Lasair%20tutorial.pdf):

1. Create a Lasair account at https://lasair.lsst.ac.uk/register/ (if you
   don't have one yet).
2. *(Optional)* Try the [Lasair quickstart](https://lasair-lsst.readthedocs.io/en/main/quickstart.html)
   to get familiar with filters and the alert format.
3. Create a **Watchmap** in Lasair, uploading the MOC FITS file you
   generated in step 3 (`output/event_visualisation/<event_id>_roi.fits`).
4. Create a **Filter** associated to that watchmap (Filter Builder →
   "Watchmaps Matching" → select your watchmap). Optionally enable email
   notifications, or just check matches later via "Run Filter".

This is what actually starts collecting alerts for your event's sky region
overnight ; make sure it's set up before you leave for the day.

## Day 2

### 5. Analyse the alerts received

Come back to Lasair and check your filter/watchmap for any matches
collected overnight:

- In Lasair, open your filter and click **Run Filter** (or check your email
  if you enabled notifications) to see the alerts matched against your
  watchmap.
- Inspect the matched objects: positions, magnitudes/fluxes, classification
  fields exposed by the Lasair schema browser.
- Cross-check candidates against the event's properties from step 3
  (distance, sky position, classification) to judge whether any are
  plausible electromagnetic counterparts versus unrelated transients
  (supernovae, variable stars, asteroids...).
