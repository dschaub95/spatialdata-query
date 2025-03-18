##
import time

import shapely
import geopandas
from spatialdata.datasets import blobs
import numpy as np
from dask.dataframe import DataFrame as DaskDataFrame
import dask.dataframe as dd
from geopandas import GeoDataFrame
from dask_geopandas import GeoDataFrame as DaskGeoDataFrame
from geoarrow.rust.core import to_geopandas
import pandas as pd
from geoarrow.rust.io import read_parquet

RNG = np.random.default_rng(2021)

def generate_blobs_sdata(n_points):
    sdata = blobs(length=100, n_points=n_points)
    return sdata

def _generate_random_points_xy(n: int, bbox: tuple[int, int]) -> np.ndarray:
    minx = bbox[0]
    maxx = bbox[1]

    return RNG.random((n, 2)) * (maxx - minx) + minx


if __name__ == "__main__":
    start = time.time()
    ga = read_parquet("random_points.parquet")
    print(f"geoarrow read parquet: {time.time() - start}")

    start = time.time()
    gdf_ga = to_geopandas(ga)
    print(f"geoarrow to geopandas: {time.time() - start}")
