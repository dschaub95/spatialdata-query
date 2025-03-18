import spatialdata as sd


def run_sdata_polygon_query(sdata, polygon):
    queried = sd.polygon_query(sdata["blobs_points"], polygon, "global")
    return queried.compute()
