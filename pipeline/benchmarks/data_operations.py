import sys
import time
import os
import pandas as pd

# Fix dask-expr issue
import dask
dask.config.set({'dataframe.query-planning': False})
import dask.dataframe as dd

# Add project root to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../.."))
sys.path.insert(0, project_root)

from pipeline.utils import save_results
from pipeline.data_generation import generate_points


# TODO: write numpy
# TODO: read numpy
# TODO: write geoparquet (geopandas)
# TODO: read geoparquet (geopanas)
# TODO: read parquet + compute (dask pandas)
# TODO: read geoparquet + compute (dask geopandas)
# TODO: to_parquet() dask dataframe
# TODO: to_geoparquet() dask geopandas
# TODO: geoarrow read_parquet (and compute it)
# TODO: geoarrow write parquet (and compute it

def benchmark_write_parquet(df, filename):
    start = time.time()
    df.to_parquet(filename)
    end = time.time()
    return end - start

def benchmark_read_parquet(filename):
    start = time.time()
    data = pd.read_parquet(filename)
    end = time.time()
    return end - start

def benchmark_dask_compute(ddf):
    start = time.time()
    ddf.compute()
    end = time.time()
    return end - start

def run_benchmark():
    n_points = int(snakemake.params.n_points)
    n_repeats = int(snakemake.params.n_repeats)
    output_file = snakemake.output[0]
    
    temp_dir = "temp_parquet"
    os.makedirs(temp_dir, exist_ok=True)
    temp_file = os.path.join(temp_dir, f"temp_{n_points}.parquet")
    
    results = []
    
    # Generate data
    print(f"Generating {n_points} points for data operations benchmark...")
    points = generate_points(n_points)
    df = pd.DataFrame(points, columns=['x', 'y'])
    ddf = dd.from_pandas(df, npartitions=max(1, min(10, n_points // 1000000)))
    
    # Run benchmarks
    benchmark_functions = {
        "write_parquet": lambda: benchmark_write_parquet(df, temp_file),
        "read_parquet": lambda: benchmark_read_parquet(temp_file),
        "dask_compute": lambda: benchmark_dask_compute(ddf)
    }
    
    for method_name, benchmark_func in benchmark_functions.items():
        print(f"Running benchmark {method_name} with {n_points} points...")
        
        for i in range(n_repeats):
            try:
                execution_time = benchmark_func()
                
                result_row = {
                    "method": method_name,
                    "n_points": n_points,
                    "execution_time": execution_time,
                    "repetition": i,
                }
                results.append(result_row)
            except Exception as e:
                print(f"Error in {method_name} benchmark: {str(e)}")
    
    # Clean up temporary files
    if os.path.exists(temp_file):
        os.remove(temp_file)
    
    # Save results
    save_results(results, output_file)

if __name__ == "__main__":
    run_benchmark()