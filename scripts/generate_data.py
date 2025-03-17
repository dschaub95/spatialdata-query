from spatialdata.datasets import blobs


def generate_blobs_sdata(n_points):
    sdata = blobs(length=100, n_points=n_points)
    return sdata
