"""Benchmarks for point queries over an axis-aligned bounding-box region.

Compares three methods against the same rectangular area:
  - run_sdata_bbox_query_points:    uses spatialdata.bounding_box_query
  - run_sdata_polygon_query_points: uses spatialdata.polygon_query with a rectangle
  - run_mpl_path:                   uses matplotlib.path with the same rectangle

Dataset selection (set before running):
    SDATA_BENCHMARK_DATASET=synthetic  (default) - blobs at multiple sizes
    SDATA_BENCHMARK_DATASET=real       - Xenium 2.0 zarr at a single size

    The real dataset can be obtained from the spatialdata docs, or by running
    download.py and then to_zarr.py in spatialdata-sandbox.
    Override the path with SDATA_REAL_DATA_PATH.

Running with pixi:
    pixi run bench-bbox-points

Running (quick, uses current Python env):
    asv run --python=same -b BoundingBoxQueryPointsBenchmark --quick --show-stderr -v

    # with real data:
    SDATA_BENCHMARK_DATASET=real asv run --python=same -b BoundingBoxQueryPointsBenchmark --quick --show-stderr -v

Running (full, multiple repeats):
    asv run --python=same -b BoundingBoxQueryPointsBenchmark --show-stderr -v

Comparing two commits:
    asv continuous --python=same main <branch> -b BoundingBoxQueryPointsBenchmark --show-stderr -v
    asv compare main <branch>
"""

from shapely.geometry import box

from utils.dataset import DATASET


class BoundingBoxQueryPointsBenchmark:
    """Benchmark point queries over the same rectangular region."""

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
        from utils.dataset import load_dataset, make_bbox

        self.sdata, self.element_key = load_dataset(dataset, "points")
        self.min_coordinate, self.max_coordinate = make_bbox(self.sdata)
        self.polygon = box(
            self.min_coordinate[0],
            self.min_coordinate[1],
            self.max_coordinate[0],
            self.max_coordinate[1],
        )

    def time_mpl_path(self, dataset) -> None:
        """Walltime for matplotlib-path rectangular polygon query."""
        from scripts.queries.mpl_path import run_mpl_path

        run_mpl_path(self.sdata, self.polygon, self.element_key)

    def time_spatialdata_bbox_query(self, dataset) -> None:
        """Walltime for spatialdata.bounding_box_query."""
        from scripts.queries.sdata_bbox_query_points import run_sdata_bbox_query_points

        run_sdata_bbox_query_points(
            self.sdata,
            self.min_coordinate,
            self.max_coordinate,
            self.element_key,
        )

    def time_spatialdata_polygon_query(self, dataset) -> None:
        """Walltime for spatialdata.polygon_query with a rectangle."""
        from scripts.queries.sdata_polygon_query_points import run_sdata_polygon_query_points

        run_sdata_polygon_query_points(self.sdata, self.polygon, self.element_key)

    def peakmem_mpl_path(self, dataset) -> None:
        """Peak memory for matplotlib-path rectangular polygon query."""
        from scripts.queries.mpl_path import run_mpl_path

        run_mpl_path(self.sdata, self.polygon, self.element_key)

    def peakmem_spatialdata_bbox_query(self, dataset) -> None:
        """Peak memory for spatialdata.bounding_box_query."""
        from scripts.queries.sdata_bbox_query_points import run_sdata_bbox_query_points

        run_sdata_bbox_query_points(
            self.sdata,
            self.min_coordinate,
            self.max_coordinate,
            self.element_key,
        )

    def peakmem_spatialdata_polygon_query(self, dataset) -> None:
        """Peak memory for spatialdata.polygon_query with a rectangle."""
        from scripts.queries.sdata_polygon_query_points import run_sdata_polygon_query_points

        run_sdata_polygon_query_points(self.sdata, self.polygon, self.element_key)
