from scripts.generate_data import generate_blobs_sdata
from scripts.mpl_path import run_mpl_path
from scripts.sdata_polygon_query import run_sdata_polygon_query
import pandas as pd
import seaborn as sns
import os
import time
import matplotlib.pyplot as plt


n_points_values = [
    int(1e5),
    int(1e6),
    int(1e7),
    # int(1e8),
]

n_repeats = 2

if __name__ == "__main__":
    methods = [run_sdata_polygon_query, run_mpl_path]
    execution_times = []

    for n_points in n_points_values:
        sdata = generate_blobs_sdata(n_points)
        polygon = sdata["blobs_polygons"].geometry.iloc[2]
        for method in methods:
            method_name = method.__name__
            for i in range(n_repeats):
                start_time = time.time()
                result = method(sdata, polygon)
                end_time = time.time()
                execution_time = end_time - start_time
                row = {
                    "method": method_name,
                    "n_points": n_points,
                    "execution_time": execution_time,
                    "repetition": i,
                }
                execution_times.append(row)

    # convert to dataframe
    execution_times_df = pd.DataFrame(execution_times)

    # plot
    sns.set_theme(style="ticks", context="paper")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        data=execution_times_df, x="n_points", y="execution_time", hue="method", ax=ax
    )
    ax.set_title("Polygon Query Benchmark")
    ax.set_xlabel("Number of Points")
    ax.set_ylabel("Execution Time (s)")
    save_dir = "figures"
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(
        os.path.join(save_dir, "benchmark_polygon_query.png"),
        dpi=300,
        bbox_inches="tight",
    )
