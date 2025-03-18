import numpy as np
from matplotlib.path import Path
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
import dask.array as da
from dask.dataframe import DataFrame as DaskDataFrame


def path_to_polygon(path):
    vertices = path.vertices
    codes = path.codes
    polygons = []
    current = []
    for v, code in zip(vertices, codes):
        if code == Path.MOVETO:
            if current:
                sub = np.array(current)
                if not np.allclose(sub[0], sub[-1]):
                    sub = np.vstack([sub, sub[0]])
                poly = Polygon(sub)
                if poly.is_valid and poly.area > 0:
                    polygons.append(poly)
            current = [v]
        else:
            current.append(v)
    if current:
        sub = np.array(current)
        if not np.allclose(sub[0], sub[-1]):
            sub = np.vstack([sub, sub[0]])
        poly = Polygon(sub)
        if poly.is_valid and poly.area > 0:
            polygons.append(poly)

    merged_polygons = merge_polygons_with_holes(polygons)

    return merged_polygons


def merge_polygons_with_holes(polygons):
    # Separate holes from outer boundaries.
    holes = []
    outer_polys = []
    for poly in polygons:
        if any(poly.within(other) for other in polygons if other != poly):
            holes.append(poly)
        else:
            outer_polys.append(poly)

    # Merge overlapping outer polygons.
    merged_outer = unary_union(outer_polys)
    if merged_outer.geom_type == "Polygon":
        merged_outer = [merged_outer]
    elif merged_outer.geom_type == "MultiPolygon":
        merged_outer = list(merged_outer)

    final_polys = []
    for outer in merged_outer:
        # For each outer, collect holes fully contained within.
        inner_holes = [hole.exterior.coords[:] for hole in holes if hole.within(outer)]
        final_polys.append(Polygon(outer.exterior.coords, inner_holes))

    return final_polys[0] if len(final_polys) == 1 else MultiPolygon(final_polys)


def polygon_to_path(polygon: Polygon | MultiPolygon) -> Path:
    if isinstance(polygon, MultiPolygon):
        raise ValueError("MultiPolygon not supported")
    elif isinstance(polygon, Polygon):
        return Path.make_compound_path(
            Path(np.asarray(polygon.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors],
        )
    else:
        raise ValueError("Invalid polygon type")


def run_mpl_path(sdata, polygon):
    return polygon_query_mpl(sdata["blobs_points"], polygon).compute()


def polygon_query_mpl(points: DaskDataFrame, polygon):
    point_coords = points[["x", "y"]].compute().values
    outer_path = Path(np.asarray(polygon.exterior.coords)[:, :2])
    inner_paths = [Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors]
    # first filter with outer path
    mask = outer_path.contains_points(point_coords)
    if len(inner_paths) > 0:
        sub_point_coords = point_coords[mask]
        # then filter with inner paths to handle holes
        sub_masks = [
            inner_path.contains_points(sub_point_coords) for inner_path in inner_paths
        ]
        # combine masks
        mask[mask] = np.logical_not(np.logical_or.reduce(sub_masks))
    # convert mask to dask array
    partition_lengths = points.map_partitions(len).compute()
    mask_dask = da.from_array(mask, chunks=partition_lengths.tolist())
    points_filtered = points.loc[mask_dask]
    return points_filtered
