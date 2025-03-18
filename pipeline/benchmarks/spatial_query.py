import sys
import time
import os
import numpy as np

# Fix dask-expr issue
import dask

dask.config.set({"dataframe.query-planning": False})

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../.."))
sys.path.insert(0, project_root)

# Create data generation functions directly in this file to avoid import issues
RNG = np.random.default_rng(2021)

# Import local utils module
from pipeline.utils import save_results
from pipeline.data_generation import generate_points, generate_polygons


def run_simple_polygon_query(points, polygon):
    """Simple polygon containment check implementation"""
    # Convert points to a numpy array if not already
    if not isinstance(points, np.ndarray):
        points = np.array(points)

    # Check polygon containment
    results = []
    for point in points:
        contains = point_in_polygon(point, polygon)
        if contains:
            results.append(point)

    return np.array(results)


def point_in_polygon(point, polygon):
    """Ray casting algorithm to determine if a point is inside a polygon"""
    x, y = point
    n = len(polygon)
    inside = False

    p1x, p1y = polygon[0]
    for i in range(1, n):
        p2x, p2y = polygon[i]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def run_spatialdata_query(points, polygon):
    """Numpy-based polygon query implementation"""
    # TODO:
    return


def run_benchmark():
    n_points = int(snakemake.params.n_points)
    n_repeats = int(snakemake.params.n_repeats)
    output_file = snakemake.output[0]

    methods = [run_simple_polygon_query, run_spatialdata_query]
    results = []

    # Generate data
    print(f"Generating {n_points} points for polygon query benchmark...")
    points = generate_points(n_points)
    polygon = generate_polygons(1)[0]  # Generate a single polygon

    # Run benchmarks
    for method in methods:
        method_name = method.__name__
        print(f"Running benchmark with method {method_name} and {n_points} points...")

        for i in range(n_repeats):
            start_time = time.time()
            _ = method(points, polygon)
            end_time = time.time()
            execution_time = end_time - start_time

            result_row = {
                "method": method_name,
                "n_points": n_points,
                "execution_time": execution_time,
                "repetition": i,
            }
            results.append(result_row)

    # Save results
    save_results(results, output_file)


if __name__ == "__main__":
    from dataclasses import dataclass, field


    @dataclass
    class Params:
        n_points: int
        n_repeats: int = 2


    @dataclass
    class SnakemakeConfig:
        params: Params
        output: list = field(init=False)

        def __post_init__(self):
            self.output = [f"results/spatial_query_results_{self.params.n_points}.csv"]


    snakemake = SnakemakeConfig(Params(10000))
    run_benchmark()
