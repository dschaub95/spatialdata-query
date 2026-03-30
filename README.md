# SpatialData Query Benchmarks

This repository contains benchmarks for spatial query operations using SpatialData
and other backends. Two toolchains are available:

- **[pixi](https://pixi.sh)** — recommended; installs the environment and exposes
  all workflows as named tasks (`pixi run <task>`)
- **[uv](https://github.com/astral-sh/uv)** — alternative; manages the virtualenv
  manually

---

## Setup

### With pixi (recommended)

```bash
pixi install
```

### With uv

```bash
uv venv && source .venv/bin/activate && uv pip install -e .
# Also install profiling tools:
uv pip install -e ".[profiling]"
```

---

## Pixi tasks

### ASV benchmarks

| Task | Saves results | What it does |
|------|:---:|-------------|
| `pixi run bench-quick` | No | Sanity check — runs every benchmark once to confirm nothing crashes |
| `pixi run bench-spatial` | Yes | Full run of `SpatialQueryBenchmark` only (3 repeats) |
| `pixi run bench-data` | Yes | Full run of `DataOperationsBenchmark` only (3 repeats) |
| `pixi run bench-full` | Yes | Full run of all benchmark suites (3 repeats) |
| `pixi run bench-report` | — | Publish saved results and open an interactive HTML report |

> **Note:** `bench-quick` does **not** save results to `.asv/results/`, so
> running it before `bench-report` will produce an empty report with 404 errors
> for all graphs. Use `bench-spatial`, `bench-data`, or `bench-full` to
> generate results that `bench-report` can display.

To tweak parameters (number of repeats, specific commits, branch comparison)
run the underlying `asv` command directly, e.g.:

```bash
# Compare HEAD against main
pixi run asv continuous --python=same main HEAD --show-stderr -v
# View comparison table
pixi run asv compare main HEAD
```

### CPU profiling with py-spy

Each query strategy has its own standalone script under `scripts/`:

| Script | Strategy |
|--------|----------|
| `scripts/profile_spatialdata_query.py` | `spatialdata.polygon_query` (default) |
| `scripts/profile_mpl_path.py` | matplotlib-path approach |

```bash
# Default — profiles the spatialdata strategy
pixi run profile-pyspy

# Profile the matplotlib-path strategy instead
PROFILE_SCRIPT=scripts/profile_mpl_path.py pixi run profile-pyspy
```

This writes `profile.speedscope.json`. Open it with:

```bash
speedscope profile.speedscope.json   # requires: npm install -g speedscope
```

> **macOS note:** py-spy may require `sudo` to attach to the process:
> ```bash
> sudo PROFILE_SCRIPT=scripts/profile_spatialdata_query.py pixi run profile-pyspy
> ```
> If `sudo` changes PATH, run py-spy directly:
> ```bash
> sudo py-spy record --gil -o profile.speedscope.json --format speedscope \
>     -- .pixi/envs/default/bin/python scripts/profile_spatialdata_query.py
> ```

To tweak the sampling rate or other py-spy flags, run it directly:

```bash
py-spy record --rate 200 --gil -o profile.speedscope.json --format speedscope \
    -- python scripts/profile_spatialdata_query.py
```

### Memory profiling with memray

```bash
# Default — profiles the spatialdata strategy; two steps: record then report
pixi run profile-memray
pixi run profile-memray-report   # opens memray-flamegraph-memray-output.html

# Profile the matplotlib-path strategy instead
PROFILE_SCRIPT=scripts/profile_mpl_path.py pixi run profile-memray
pixi run profile-memray-report
```

To tweak the report format or output path, run memray directly:

```bash
memray run -o memray-output.bin scripts/profile_spatialdata_query.py
memray flamegraph memray-output.bin              # non-temporal flamegraph
memray flamegraph --temporal memray-output.bin   # temporal flamegraph
memray summary memray-output.bin                 # text summary
```

---

## Benchmarks

| File | Class | What it measures |
|------|-------|-----------------|
| `benchmarks/benchmark_spatial_query.py` | `SpatialQueryBenchmark` | Walltime and peak memory of `spatialdata.polygon_query` vs. the matplotlib-path method at 1e5 and 1e6 points |
| `benchmarks/benchmark_data_operations.py` | `DataOperationsBenchmark` | Walltime and peak memory for pandas `to_parquet`, `read_parquet`, and dask `compute` at 1e3–1e6 points |

Each class follows the ASV convention:
- `time_*` — wall-clock timing
- `peakmem_*` — peak RSS memory
- `params` / `param_names` — parameterised over `n_points`
- `setup` / `teardown` — per-run data generation and cleanup (not timed)

### Profiling scripts

`scripts/benchmark_polygon_query.py` is a standalone script used by the
py-spy and memray tasks. It benchmarks the same two polygon-query methods
(`run_sdata_polygon_query` and `run_mpl_path`) and verifies that both return
identical results. Edit `n_points_values` and `n_repeats` at the top of the
file to adjust the dataset size and repetition count before profiling.
