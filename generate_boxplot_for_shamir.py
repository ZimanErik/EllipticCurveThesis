import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 generate_boxplot_for_shamir.py <file.csv>")
        sys.exit(1)

    file_path = sys.argv[1]
    if not os.path.isfile(file_path):
        print(f"File not found: {file_path}")
        sys.exit(1)

    try:
        df_raw = pd.read_csv(file_path, header=None, sep=r',+', engine='python')
    except pd.errors.EmptyDataError:
        print(f"Error: CSV file '{file_path}' is empty or malformed")
        sys.exit(1)
    if df_raw.iloc[:, -1].isna().all():
        df_raw = df_raw.iloc[:, :-1]

    if df_raw.shape[0] != 2 or df_raw.shape[1] != 1000:
        print(f"CSV must contain exactly 2 rows and 1000 columns. Found {df_raw.shape[0]} rows and {df_raw.shape[1]} columns.")
        sys.exit(1)

    labels = ["Shamir's trick", "Naive scalar multiplication"]
    df_melted = pd.DataFrame({
        "value": df_raw.iloc[0].tolist() + df_raw.iloc[1].tolist(),
        " ": [labels[0]] * 1000 + [labels[1]] * 1000
    })

    for i, label in enumerate(labels):
        group = df_melted[df_melted[" "] == label]["value"]
        min_val = np.min(group)
        q1 = np.percentile(group, 25)
        median = np.median(group)
        mean = np.mean(group)
        q3 = np.percentile(group, 75)
        max_val = np.max(group)
        std = np.std(group)

        print(f"Statistics for {label}:")
        print(f"   Min: {min_val:.0f}")
        print(f"   1st Quartile (Q1): {q1:.0f}")
        print(f"   Median: {median:.0f}")
        print(f"   Mean: {mean:.0f}")
        print(f"   3rd Quartile (Q3): {q3:.0f}")
        print(f"   Max: {max_val:.0f}")
        print(f"   Standard Deviation: {std:.0f}")
        print("-" * 40)

    plt.figure(figsize=(8, 6))
    ax = sns.boxplot(
        x=" ", y="value", data=df_melted,
        flierprops=dict(marker='o', markersize=2, alpha=0.2)
    )

    plt.title("Shamir's Trick vs Naive Scalar Point Multiplication (for 256-bit curves)")
    plt.ylabel("Execution Time (in μs)")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
