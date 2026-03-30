import numpy as np
import pandas as pd


def write_parquet(df: pd.DataFrame, path: str) -> None:
    df.to_parquet(path)


def read_parquet(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


def dask_compute(ddf):
    return ddf.compute()


if __name__ == "__main__":
    import tempfile
    import time
    from pathlib import Path

    import dask.dataframe as dd

    n_points_values = [
        int(1e3),
        int(1e4),
        int(1e5),
        int(1e6),
    ]
    n_repeats = 2

    with tempfile.TemporaryDirectory() as tmpdir:
        for n_points in n_points_values:
            rng = np.random.default_rng(2021)
            points = rng.random((n_points, 2)) * 100
            df = pd.DataFrame(points, columns=["x", "y"])
            ddf = dd.from_pandas(df, npartitions=max(1, min(10, n_points // 1_000_000)))
            parquet_path = str(Path(tmpdir) / "data.parquet")
            df.to_parquet(parquet_path)

            for op_name, fn in [
                ("write_parquet", lambda: write_parquet(df, str(Path(tmpdir) / "bench.parquet"))),
                ("read_parquet", lambda: read_parquet(parquet_path)),
                ("dask_compute", lambda: dask_compute(ddf)),
            ]:
                for i in range(n_repeats):
                    t0 = time.perf_counter()
                    fn()
                    elapsed = time.perf_counter() - t0
                    print(f"n_points={n_points:>8}  op={op_name}  rep={i}  time={elapsed:.4f}s")
