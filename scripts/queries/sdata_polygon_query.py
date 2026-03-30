import spatialdata as sd


def run_sdata_polygon_query(sdata, polygon):
    queried = sd.polygon_query(sdata["blobs_points"], polygon, "global")
    return queried.compute()


if __name__ == "__main__":
    import time

    from utils.generate_data import add_circle_polygon, generate_blobs_sdata

    n_points_values = [
        int(1e5),
        int(1e6),
        # int(1e7),
    ]
    n_repeats = 2

    for n_points in n_points_values:
        sdata = generate_blobs_sdata(n_points)
        polygon = add_circle_polygon(sdata, has_hole=True, radius=10.1)
        for i in range(n_repeats):
            t0 = time.perf_counter()
            result = run_sdata_polygon_query(sdata, polygon)
            elapsed = time.perf_counter() - t0
            print(f"n_points={n_points:>8}  rep={i}  time={elapsed:.3f}s  result={result.shape}")
