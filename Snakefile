# type: ignore

import os

# Configuration
N_REPEATS = 3
N_POINTS_VALUES = [
    int(1e3),
    int(1e4),
    int(1e5),
    int(1e6), 
    # int(1e7),
    # int(1e8),
]

# Create required directories
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)
os.makedirs("pipeline/benchmarks", exist_ok=True)

# Define final output files
rule all:
    input:
        "figures/spatial_query_benchmark.png",
        "figures/data_operations_benchmark.png"

# Run the spatial query benchmark with multiple repetitions
rule run_spatial_query_benchmark:
    output:
        "results/spatial_query_results_{n_points}.csv"
    params:
        n_points = lambda wildcards: wildcards.n_points,
        n_repeats = N_REPEATS
    script:
        "pipeline/benchmarks/spatial_query.py"

# Run the parquet operations benchmark with multiple repetitions
rule run_data_operations_benchmark:
    output:
        "results/data_operations_results_{n_points}.csv"
    params:
        n_points = lambda wildcards: wildcards.n_points,
        n_repeats = N_REPEATS
    script:
        "pipeline/benchmarks/data_operations.py"

# Collect and plot the polygon query benchmark results
rule plot_spatial_query:
    input:
        expand("results/spatial_query_results_{n_points}.csv", n_points=N_POINTS_VALUES)
    output:
        "figures/spatial_query_benchmark.png"
    script:
        "pipeline/benchmarks/plot_spatial_query.py"

# Collect and plot the data operations benchmark results
rule plot_data_operations:
    input:
        expand("results/data_operations_results_{n_points}.csv",
            n_points=N_POINTS_VALUES)
    output:
        "figures/data_operations_benchmark.png"
    script:
        "pipeline/benchmarks/plot_data_operations.py"