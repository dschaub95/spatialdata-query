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
PROFILE_SCRIPT=scripts/profile_mpl_path.py pixi run profile-pyspy   # mpl-path strategy

pixi run profile-memray       # memory recording → memray-output.bin  (spatialdata strategy)
PROFILE_SCRIPT=scripts/profile_mpl_path.py pixi run profile-memray  # mpl-path strategy
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

## Tasks

### Benchmarks (ASV)

| Task | Saves results | Description |
|------|:---:|-------------|
| `bench-quick` | No | Crash-check — one run, no results written |
| `bench-spatial` | Yes | `SpatialQueryBenchmark` — 3 repeats |
| `bench-data` | Yes | `DataOperationsBenchmark` — 3 repeats |
| `bench-full` | Yes | All suites — 3 repeats |
| `bench-report` | — | Print latest results to terminal |

Branch comparison:
```bash
pixi run asv continuous --python=same main HEAD --show-stderr -v
pixi run asv compare main HEAD
```

### Profiling

Edit `n_points_values` and `n_repeats` at the top of the profile scripts to tune the workload.
macOS py-spy may require `sudo`. Speedscope requires `npm install -g speedscope`.
