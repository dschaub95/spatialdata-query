import numpy as np
from matplotlib.path import Path


def run_mpl_path(sdata, polygon):
    coords = sdata["blobs_points"].compute()[["x", "y"]].values
    path = Path.make_compound_path(
        Path(np.asarray(polygon.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors],
    )
    mask = path.contains_points(coords)
    points_filtered = sdata["blobs_points"].compute().loc[mask]
    return points_filtered
