"""Benchmarks for spatial polygon query operations.

Compares two methods for querying points within a polygon:
  - run_sdata_polygon_query: uses spatialdata.polygon_query
  - run_mpl_path:            uses matplotlib.path for point-in-polygon

Running (quick, uses current Python env):
    asv run --python=same -b SpatialQueryBenchmark --quick --show-stderr -v

Running (full, multiple repeats):
    asv run --python=same -b SpatialQueryBenchmark --show-stderr -v

Comparing two commits:
    asv continuous --python=same main <branch> -b SpatialQueryBenchmark --show-stderr -v
    asv compare main <branch>

HTML report:
    asv publish && asv preview
"""


class SpatialQueryBenchmark:
    """Benchmark spatial polygon query methods at different point counts."""

    timeout = 600
    repeat = 3
    number = 1
    warmup_time = 0
    processes = 1

    params = [[int(1e5), int(1e6)]]
    param_names = ["n_points"]

    def setup(self, n_points: int) -> None:
        from scripts.generate_data import add_circle_polygon, generate_blobs_sdata

        self.sdata = generate_blobs_sdata(n_points)
        self.polygon = add_circle_polygon(self.sdata, has_hole=True, radius=10.1)

    def time_mpl_path(self, n_points: int) -> None:
        """Walltime for matplotlib-path polygon query."""
        from scripts.mpl_path import run_mpl_path

        run_mpl_path(self.sdata, self.polygon)

    def time_spatialdata_query(self, n_points: int) -> None:
        """Walltime for spatialdata.polygon_query."""
        from scripts.sdata_polygon_query import run_sdata_polygon_query

        run_sdata_polygon_query(self.sdata, self.polygon)

    def peakmem_mpl_path(self, n_points: int) -> None:
        """Peak memory for matplotlib-path polygon query."""
        from scripts.mpl_path import run_mpl_path

        run_mpl_path(self.sdata, self.polygon)

    def peakmem_spatialdata_query(self, n_points: int) -> None:
        """Peak memory for spatialdata.polygon_query."""
        from scripts.sdata_polygon_query import run_sdata_polygon_query

        run_sdata_polygon_query(self.sdata, self.polygon)
