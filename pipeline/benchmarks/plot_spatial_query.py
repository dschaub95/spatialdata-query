import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# Fix dask-expr issue
import dask
dask.config.set({'dataframe.query-planning': False})

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../.."))
sys.path.insert(0, project_root)

from pipeline.utils import aggregate_results, setup_sns_style

def plot_polygon_query_results():
    # Concatenate all results
    results_dfs = []
    for input_file in snakemake.input:
        df = pd.read_csv(input_file)
        results_dfs.append(df)
    
    all_results = pd.concat(results_dfs)
    
    # Aggregate results to get mean and std
    aggregated = aggregate_results(all_results)
    
    # Create the plot
    setup_sns_style()
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot mean with error bars showing standard deviation
    for method in aggregated['method'].unique():
        method_data = aggregated[aggregated['method'] == method]
        ax.errorbar(
            x=method_data['n_points'], 
            y=method_data['execution_time_mean'],
            yerr=method_data['execution_time_std'],
            label=method,
            marker='o',
            capsize=4
        )
    
    ax.set_title("Polygon Query Benchmark")
    ax.set_xlabel("Number of Points")
    ax.set_ylabel("Execution Time (s)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(snakemake.output[0], dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_polygon_query_results()