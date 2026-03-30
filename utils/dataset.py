"""Dataset configuration and loading for benchmarks and profiling scripts.

Controls which dataset is used via environment variables:

    SDATA_BENCHMARK_DATASET=synthetic  (default) — generate blobs at a given n_points
    SDATA_BENCHMARK_DATASET=real       — load the Xenium 2.0 zarr from SDATA_REAL_DATA_PATH

    The real dataset can be obtained from the spatialdata docs, or by running
    download.py and then to_zarr.py in spatialdata-sandbox.

ASV orchestrates the different synthetic sizes; profiling scripts use a single fixed size.

Element types
-------------
element_type="points"  — real: "transcripts",      synthetic: "blobs_points"
element_type="shapes"  — real: "cell_boundaries",  synthetic: "blobs_polygons"
"""

import os

DATASET = os.environ.get("SDATA_BENCHMARK_DATASET", "synthetic")
REAL_DATA_PATH = os.environ.get(
    "SDATA_REAL_DATA_PATH",
    "/Users/macbook/embl/projects/basel/spatialdata-sandbox/xenium_2.0.0_io/data.zarr",
)

_ELEMENT_KEYS = {
    "real":      {"points": "transcripts",   "shapes": "cell_boundaries"},
    "synthetic": {"points": "blobs_points",  "shapes": "blobs_polygons"},
}


def load_dataset(dataset, element_type):
    """Load and subset a SpatialData object to a single element.

    Parameters
    ----------
    dataset : "real" or int
        "real" loads the Xenium 2.0 zarr; an int generates blobs with that many points.
    element_type : "points" or "shapes"

    Returns
    -------
    sdata : SpatialData  — subset to the single element
    element_key : str    — key into sdata to query against
    """
    import spatialdata as sd

    from utils.generate_data import generate_blobs_sdata

    kind = "real" if dataset == "real" else "synthetic"
    element_key = _ELEMENT_KEYS[kind][element_type]

    if dataset == "real":
        sdata = sd.read_zarr(REAL_DATA_PATH)
    else:
        sdata = generate_blobs_sdata(dataset)

    sdata = sdata.subset([element_key])
    return sdata, element_key


def make_polygon(sdata):
    """Circle-with-hole polygon centred on the sdata extent (radius = 10 % of shorter axis)."""
    import spatialdata as sd

    from utils.generate_data import add_circle_polygon

    extent = sd.get_extent(sdata)
    xrange = extent["x"][1] - extent["x"][0]
    yrange = extent["y"][1] - extent["y"][0]
    radius = min(xrange, yrange) * 0.1
    return add_circle_polygon(sdata, has_hole=True, radius=radius)


def make_bbox(sdata):
    """Bounding box covering the central 20 % of the sdata extent.

    Returns
    -------
    min_coordinate : list[float]  — [x_min, y_min]
    max_coordinate : list[float]  — [x_max, y_max]
    """
    import spatialdata as sd

    extent = sd.get_extent(sdata)
    xmin, xmax = extent["x"]
    ymin, ymax = extent["y"]
    cx = (xmin + xmax) / 2
    cy = (ymin + ymax) / 2
    xrange = xmax - xmin
    yrange = ymax - ymin
    return (
        [cx - xrange * 0.1, cy - yrange * 0.1],
        [cx + xrange * 0.1, cy + yrange * 0.1],
    )
