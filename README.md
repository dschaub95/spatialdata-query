# SpatialData Query Benchmarks

This repository contains benchmarks for spatial query operations using SpatialData 
and other backends.

## Setup

1. Install dependencies:

Install dependencies using `uv`:

```bash
uv venv
source .venv/bin/activate
```

2. Run the benchmarks:

```bash
snakemake --cores all
```

## Benchmarks

The benchmarks are organized into two categories:

1. **Data handling operations**: Tests the performance of reading and writing 
   vector geometries into different file formats, or converting between different 
   in-memory representations, using different libraries.
2. **Spatial queries**: Tests the performance of spatial query operations using 
   different methods.

Each benchmark is run 3 times to calculate average performance and standard deviation.

## Outputs

The results are saved in the following locations:

- Raw benchmark results: `results/`
- Plots: `figures/`