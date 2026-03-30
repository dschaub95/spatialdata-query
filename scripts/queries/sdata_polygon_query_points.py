import spatialdata as sd


def run_sdata_polygon_query_points(sdata, polygon, element_key):
    queried = sd.polygon_query(sdata[element_key], polygon, "global")
    return queried.compute()


if __name__ == "__main__":
    import time

    from utils.dataset import DATASET, load_dataset, make_polygon

    dataset = "real" if DATASET == "real" else int(1e6)

    sdata, element_key = load_dataset(dataset, "points")
    polygon = make_polygon(sdata)
    t0 = time.perf_counter()
    result = run_sdata_polygon_query_points(sdata, polygon, element_key)
    elapsed = time.perf_counter() - t0
    print(f"dataset={dataset!r}  time={elapsed:.3f}s  result={result.shape}")
