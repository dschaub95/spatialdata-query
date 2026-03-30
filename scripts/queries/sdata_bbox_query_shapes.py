import spatialdata as sd


def run_sdata_bbox_query_shapes(sdata, min_coordinate, max_coordinate, element_key):
    return sd.bounding_box_query(
        sdata[element_key],
        axes=("x", "y"),
        min_coordinate=min_coordinate,
        max_coordinate=max_coordinate,
        target_coordinate_system="global",
    )


if __name__ == "__main__":
    import time

    from utils.dataset import DATASET, load_dataset, make_bbox

    dataset = "real" if DATASET == "real" else int(1e6)

    sdata, element_key = load_dataset(dataset, "shapes")
    min_coordinate, max_coordinate = make_bbox(sdata)
    t0 = time.perf_counter()
    result = run_sdata_bbox_query_shapes(sdata, min_coordinate, max_coordinate, element_key)
    elapsed = time.perf_counter() - t0
    print(f"dataset={dataset!r}  time={elapsed:.3f}s  result={result.shape}")
