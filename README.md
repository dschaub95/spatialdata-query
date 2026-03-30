# SpatialData Query Benchmarks

Benchmarks for spatial polygon query operations using [SpatialData](https://spatialdata.scverse.org) and alternative backends.

---

## Quick start

```bash
pixi install                  # set up environment

pixi run bench-spatial        # benchmark spatial query (saves results)
pixi run bench-data           # benchmark data I/O (saves results)
pixi run bench-report         # show latest results in terminal

pixi run profile-pyspy        # CPU flame graph  → profile.speedscope.json  (spatialdata strategy)
PROFILE_SCRIPT=scripts/queries/mpl_path.py pixi run profile-pyspy   # mpl-path strategy

pixi run profile-memray       # memory recording → memray-output.bin  (spatialdata strategy)
PROFILE_SCRIPT=scripts/queries/mpl_path.py pixi run profile-memray  # mpl-path strategy
pixi run profile-memray-report  # flame graph    → opens in browser
```

> `bench-quick` is a crash-check only — it does **not** save results.

---

## Setup

```bash
pixi install          # recommended (conda-forge + PyPI, all deps)

# or with uv:
uv venv && source .venv/bin/activate
uv pip install -e ".[profiling]"
```

---

## Dependencies

By default `spatialdata` is installed from PyPI.

To benchmark against a **local checkout** (e.g. a dev branch):

**pixi** — add it to `[tool.pixi.pypi-dependencies]` in `pyproject.toml` (pixi overrides the PyPI version), replacing the path with the actual location of your local `spatialdata` clone:

```toml
[tool.pixi.pypi-dependencies]
spatialdata = { path = "/path/to/your/spatialdata", editable = true }
```

Then re-run `pixi install`.

**uv / venv** — install the local version on top of the existing environment, replacing the path with the actual location of your local clone:

```bash
uv pip install -e /path/to/your/spatialdata
# or
pip install -e /path/to/your/spatialdata
```

To go back to the PyPI version, remove the override from `pyproject.toml` and re-run `pixi install`, or run `pip install spatialdata` in the venv.

---

## Tasks

### Benchmarks (ASV)

| Task | Saves results | Description |
|------|:---:|-------------|
| `bench-quick` | No | Crash-check — one run, no results written |
| `bench-spatial` | Yes | `SpatialQueryBenchmark` — 3 repeats |
| `bench-data` | Yes | `DataOperationsBenchmark` — 3 repeats |
| `bench-full` | Yes | All suites — 3 repeats |
| `bench-report` | — | Print latest results to terminal |

#### Dataset selection for `SpatialQueryBenchmark`

By default the benchmark generates **synthetic blobs** at multiple point counts (`1e5`, `1e6`).

Set `SDATA_BENCHMARK_DATASET=real` to benchmark against the **Xenium 2.0 zarr** instead — a single run at real dataset size (~12 M transcripts).

```bash
# synthetic (default) — multiple sizes
pixi run bench-spatial

# real data — single size
SDATA_BENCHMARK_DATASET=real pixi run bench-spatial

# override the zarr path (defaults to spatialdata-sandbox/xenium_2.0.0_io/data.zarr)
SDATA_BENCHMARK_DATASET=real SDATA_REAL_DATA_PATH=/path/to/data.zarr pixi run bench-spatial
```

The Xenium dataset can be obtained from the [spatialdata docs](https://spatialdata.scverse.org) or by running `download.py` and then `to_zarr.py` in `spatialdata-sandbox`.

Branch comparison:
```bash
pixi run asv continuous --python=same main HEAD --show-stderr -v
pixi run asv compare main HEAD
```

### Profiling

Edit `n_points_values` and `n_repeats` at the top of the profile scripts to tune the workload.
macOS py-spy may require `sudo`. Speedscope requires `npm install -g speedscope`.
