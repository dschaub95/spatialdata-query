"""Dataset configuration and loading for benchmarks and profiling scripts.

Controls which dataset is used via environment variables:

    SDATA_BENCHMARK_DATASET=synthetic  (default) — generate blobs at a given n_points
    SDATA_BENCHMARK_DATASET=real       — load the Xenium 2.0 zarr from SDATA_REAL_DATA_PATH

    The real dataset can be obtained from the spatialdata docs, or by running
    download.py and then to_zarr.py in spatialdata-sandbox.

ASV orchestrates the different synthetic sizes; profiling scripts loop over them directly.
"""

import os

DATASET = os.environ.get("SDATA_BENCHMARK_DATASET", "synthetic")
REAL_DATA_PATH = os.environ.get(
    "SDATA_REAL_DATA_PATH",
    "/Users/macbook/embl/projects/basel/spatialdata-sandbox/xenium_2.0.0_io/data.zarr",
)


def load_dataset(dataset):
    """Load sdata, polygon, and points_key for a given dataset spec.

    Parameters
    ----------
    dataset : "real" or int
        "real" loads the Xenium 2.0 zarr; an int generates blobs with that many points.

    Returns
    -------
    sdata : SpatialData
    polygon : shapely geometry  — query polygon centred on the data extent
    points_key : str            — key into sdata.points to query against
    """
    import spatialdata as sd

    from utils.generate_data import add_circle_polygon, generate_blobs_sdata

    if dataset == "real":
        sdata = sd.read_zarr(REAL_DATA_PATH)
        extent = sd.get_extent(sdata)
        xrange = extent["x"][1] - extent["x"][0]
        yrange = extent["y"][1] - extent["y"][0]
        radius = min(xrange, yrange) * 0.1
        points_key = "transcripts"
    else:
        sdata = generate_blobs_sdata(dataset)
        radius = 10.1
        points_key = "blobs_points"

    polygon = add_circle_polygon(sdata, has_hole=True, radius=radius)
    return sdata, polygon, points_key
