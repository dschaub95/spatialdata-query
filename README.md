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

| Task | What it does |
|------|-------------|
| `pixi run bench-quick` | Run all benchmarks once (quick sanity check) |
| `pixi run bench-spatial` | Quick run of `SpatialQueryBenchmark` only |
| `pixi run bench-data` | Quick run of `DataOperationsBenchmark` only |
| `pixi run bench-full` | Full run (3 repeats) — results saved to `.asv/results/` |
| `pixi run bench-report` | Publish results and open an interactive HTML report |

To tweak parameters (number of repeats, specific commits, branch comparison)
run the underlying `asv` command directly, e.g.:

```bash
# Compare HEAD against main
pixi run asv continuous --python=same main HEAD --show-stderr -v
# View comparison table
pixi run asv compare main HEAD
```

### CPU profiling with py-spy

```bash
# Default (profiles scripts/benchmark_polygon_query.py)
pixi run profile-pyspy

# Override the script via env var
PROFILE_SCRIPT=scripts/other.py pixi run profile-pyspy
```

This writes `profile.speedscope.json`. Open it with:

```bash
speedscope profile.speedscope.json   # requires: npm install -g speedscope
```

> **macOS note:** py-spy may require `sudo` to attach to the process:
> ```bash
> sudo PROFILE_SCRIPT=scripts/benchmark_polygon_query.py pixi run profile-pyspy
> ```
> If `sudo` changes PATH, run py-spy directly:
> ```bash
> sudo py-spy record --gil -o profile.speedscope.json --format speedscope \
>     -- .pixi/envs/default/bin/python scripts/benchmark_polygon_query.py
> ```

To tweak the sampling rate or other py-spy flags, run it directly:

```bash
py-spy record --rate 200 --gil -o profile.speedscope.json --format speedscope \
    -- python scripts/benchmark_polygon_query.py
```

### Memory profiling with memray

Two steps — record, then report:

```bash
# Step 1: record (default script; override with PROFILE_SCRIPT=other.py)
pixi run profile-memray

# Step 2: generate the HTML temporal flame graph
pixi run profile-memray-report
```

Output: `memray-flamegraph-memray-output.html`. Open in a browser:

```bash
open memray-flamegraph-memray-output.html   # macOS
xdg-open memray-flamegraph-memray-output.html  # Linux
```

To tweak the report format or output path, run memray directly:

```bash
memray run -o memray-output.bin scripts/benchmark_polygon_query.py
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
