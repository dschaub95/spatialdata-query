from spatialdata.datasets import blobs
import spatialdata as sd
import shapely


def generate_blobs_sdata(n_points):
    sdata = blobs(length=100, n_points=n_points)
    return sdata


def add_circle_polygon(sdata, radius=10, has_hole=False):
    extent = sd.get_extent(sdata)
    xmin, xmax = extent["x"]
    ymin, ymax = extent["y"]

    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    if has_hole:
        outer_circle = shapely.Point(center_x, center_y).buffer(radius)
        inner_circle = shapely.Point(center_x, center_y).buffer(radius / 2)
        circle = outer_circle.difference(inner_circle)
    else:
        circle = shapely.Point(center_x, center_y).buffer(radius)

    return circle
