import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def save_results(results, output_file):
    """Save benchmark results to a CSV file"""
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    return df

def aggregate_results(df):
    """Aggregate results to get mean and std for each method and n_points"""
    aggregated = df.groupby(['method', 'n_points']).agg({
        'execution_time': ['mean', 'std']
    }).reset_index()
    
    # Flatten the MultiIndex columns
    aggregated.columns = ['_'.join(col).strip('_') for col in aggregated.columns.values]
    return aggregated

def setup_sns_style():
    """Set up the seaborn plotting style"""
    sns.set_theme(style="ticks", context="paper")