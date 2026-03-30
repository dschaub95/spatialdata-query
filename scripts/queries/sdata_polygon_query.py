import spatialdata as sd


def run_sdata_polygon_query(sdata, polygon, points_key):
    queried = sd.polygon_query(sdata[points_key], polygon, "global")
    return queried.compute()


if __name__ == "__main__":
    import time

    from utils.dataset import DATASET, load_dataset

    dataset = "real" if DATASET == "real" else int(1e6)

    sdata, polygon, points_key = load_dataset(dataset)
    t0 = time.perf_counter()
    result = run_sdata_polygon_query(sdata, polygon, points_key)
    elapsed = time.perf_counter() - t0
    print(f"dataset={dataset!r}  time={elapsed:.3f}s  result={result.shape}")
