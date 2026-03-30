"""Benchmarks for spatial polygon query operations.

Compares two methods for querying points within a polygon:
  - run_sdata_polygon_query: uses spatialdata.polygon_query
  - run_mpl_path:            uses matplotlib.path for point-in-polygon

Dataset selection (set before running):
    SDATA_BENCHMARK_DATASET=synthetic  (default) — blobs at multiple sizes
    SDATA_BENCHMARK_DATASET=real       — Xenium 2.0 zarr at a single size

    The real dataset can be obtained from the spatialdata docs, or by running
    download.py and then to_zarr.py in spatialdata-sandbox.
    Override the path with SDATA_REAL_DATA_PATH.

Running (quick, uses current Python env):
    asv run --python=same -b SpatialQueryBenchmark --quick --show-stderr -v

    # with real data:
    SDATA_BENCHMARK_DATASET=real asv run --python=same -b SpatialQueryBenchmark --quick --show-stderr -v

Running (full, multiple repeats):
    asv run --python=same -b SpatialQueryBenchmark --show-stderr -v

Comparing two commits:
    asv continuous --python=same main <branch> -b SpatialQueryBenchmark --show-stderr -v
    asv compare main <branch>
"""

from utils.dataset import DATASET


class SpatialQueryBenchmark:
    """Benchmark spatial polygon query methods.

    With SDATA_BENCHMARK_DATASET=synthetic (default): parameterised over point counts.
    With SDATA_BENCHMARK_DATASET=real: single run against the Xenium 2.0 zarr dataset.
    """

    timeout = 600
    repeat = 3
    number = 1
    warmup_time = 0
    processes = 1

    if DATASET == "real":
        params = [["real"]]
    else:
        params = [[int(1e5), int(1e6)]]
    param_names = ["dataset"]

    def setup(self, dataset) -> None:
        from utils.dataset import load_dataset

        self.sdata, self.polygon, self.points_key = load_dataset(dataset)

    def time_mpl_path(self, dataset) -> None:
        """Walltime for matplotlib-path polygon query."""
        from scripts.queries.mpl_path import run_mpl_path

        run_mpl_path(self.sdata, self.polygon, self.points_key)

    def time_spatialdata_query(self, dataset) -> None:
        """Walltime for spatialdata.polygon_query."""
        from scripts.queries.sdata_polygon_query import run_sdata_polygon_query

        run_sdata_polygon_query(self.sdata, self.polygon, self.points_key)

    def peakmem_mpl_path(self, dataset) -> None:
        """Peak memory for matplotlib-path polygon query."""
        from scripts.queries.mpl_path import run_mpl_path

        run_mpl_path(self.sdata, self.polygon, self.points_key)

    def peakmem_spatialdata_query(self, dataset) -> None:
        """Peak memory for spatialdata.polygon_query."""
        from scripts.queries.sdata_polygon_query import run_sdata_polygon_query

        run_sdata_polygon_query(self.sdata, self.polygon, self.points_key)
