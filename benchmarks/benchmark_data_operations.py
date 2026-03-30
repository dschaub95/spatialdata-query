"""Benchmarks for data I/O operations.

Measures performance of reading and writing point data in different formats
(Pandas parquet, Dask parquet) at varying dataset sizes.

Running (quick, uses current Python env):
    asv run --python=same -b DataOperationsBenchmark --quick --show-stderr -v

Running (full, multiple repeats):
    asv run --python=same -b DataOperationsBenchmark --show-stderr -v

Comparing two commits:
    asv continuous --python=same main <branch> -b DataOperationsBenchmark --show-stderr -v
    asv compare main <branch>

HTML report:
    asv publish && asv preview
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


class DataOperationsBenchmark:
    """Benchmark parquet and dask I/O at different point counts."""

    timeout = 300
    repeat = 3
    number = 1
    warmup_time = 0
    processes = 1

    params = [[int(1e3), int(1e4), int(1e5), int(1e6)]]
    param_names = ["n_points"]

    def setup(self, n_points: int) -> None:
        import dask.dataframe as dd

        rng = np.random.default_rng(2021)
        points = rng.random((n_points, 2)) * 100
        self.df = pd.DataFrame(points, columns=["x", "y"])
        self.ddf = dd.from_pandas(self.df, npartitions=max(1, min(10, n_points // 1_000_000)))

        self._tmpdir = tempfile.TemporaryDirectory()
        self.parquet_path = str(Path(self._tmpdir.name) / "data.parquet")
        # Pre-write file so read benchmarks can run independently
        self.df.to_parquet(self.parquet_path)

    def teardown(self, n_points: int) -> None:
        self._tmpdir.cleanup()

    # ---- timing benchmarks ------------------------------------------------

    def time_write_parquet(self, n_points: int) -> None:
        """Walltime for pandas DataFrame.to_parquet."""
        path = str(Path(self._tmpdir.name) / "write_bench.parquet")
        self.df.to_parquet(path)

    def time_read_parquet(self, n_points: int) -> None:
        """Walltime for pandas read_parquet."""
        pd.read_parquet(self.parquet_path)

    def time_dask_compute(self, n_points: int) -> None:
        """Walltime for dask DataFrame.compute (in-memory partition concat)."""
        self.ddf.compute()

    # ---- peak-memory benchmarks -------------------------------------------

    def peakmem_write_parquet(self, n_points: int) -> None:
        """Peak memory for pandas DataFrame.to_parquet."""
        path = str(Path(self._tmpdir.name) / "write_bench_mem.parquet")
        self.df.to_parquet(path)

    def peakmem_read_parquet(self, n_points: int) -> None:
        """Peak memory for pandas read_parquet."""
        pd.read_parquet(self.parquet_path)

    def peakmem_dask_compute(self, n_points: int) -> None:
        """Peak memory for dask DataFrame.compute."""
        self.ddf.compute()
