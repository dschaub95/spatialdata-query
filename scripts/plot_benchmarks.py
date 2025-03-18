import os
import pandas as pd
import matplotlib.pyplot as plt

def gather_data():
    data = []
    for file in os.listdir('data'):
        if file.endswith('.csv'):
            df = pd.read_csv(os.path.join('data', file))
            data.append(df)
    return pd.concat(data, ignore_index=True)

def plot_benchmarks(data):
    for benchmark in data['benchmark'].unique():
        plt.figure()
        for data_type in data['data_type'].unique():
            subset = data[(data['benchmark'] == benchmark) & (data['data_type'] == data_type)]
            plt.plot(subset['size'], subset['duration'], label=data_type)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Size')
        plt.ylabel('Duration (s)')
        plt.title(f'Benchmark: {benchmark}')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'plots/{benchmark}.png')

if __name__ == "__main__":
    os.makedirs('plots', exist_ok=True)
    data = gather_data()
    plot_benchmarks(data)