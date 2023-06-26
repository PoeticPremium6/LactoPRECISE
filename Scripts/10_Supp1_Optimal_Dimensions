# Filter the expression matrix with the correlation of samples within a condition group.
# Average all the samples if one of the pairs in group does not pass the correlation threshold

import os
import pandas as pd
import numpy as np
import argparse
from itertools import combinations
from scipy.stats import pearsonr

def main(input_file, group_file):
    # Read the input expression matrix file
    df = pd.read_csv(input_file, index_col=0)

    # Read the group file
    groups = pd.read_csv(group_file)

    # Group columns by the same "str" before Rep_1, Rep_2, and Rep_3
    groups_dict = {}
    for _, row in groups.iterrows():
        group_name = row[1].rsplit('_', 2)[0]
        column_name = row[0]
        if group_name not in groups_dict:
            groups_dict[group_name] = []
        groups_dict[group_name].append(column_name)

    # Create output directory
    output_dir = f"Corr_Average_Filtered_{os.path.splitext(os.path.basename(input_file))[0]}"
    os.makedirs(output_dir, exist_ok=True)

    # Filter the original input expression matrix based on the thresholds
    for threshold in [0.8, 0.85, 0.9, 0.95]:
        filtered_df = df.copy()
        for group, columns in groups_dict.items():
            should_average = False
            for col1, col2 in combinations(columns, 2):
                if col1 not in df.columns or col2 not in df.columns:
                    continue
                if col1 == df.columns[0] or col1 == df.columns[1] or col2 == df.columns[0] or col2 == df.columns[1]:
                    continue
                x = df.loc[:, col1].values
                y = df.loc[:, col2].values
                r, _ = pearsonr(x, y)
                if r < threshold:
                    should_average = True
                    break
            if should_average:
                average_column = df[columns].mean(axis=1)
                for col in columns[1:]:
                    filtered_df.drop(col, axis=1, inplace=True)
                filtered_df[columns[0]] = average_column
                filtered_df.rename(columns={columns[0]: f"{columns[0]}_Average"}, inplace=True)
        output_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(input_file))[0]}_filtered_{threshold}.csv")
        filtered_df.to_csv(output_file, index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Average all the columns in a group if one pair of columns within the group does not pass the threshold and store the result in the first column of that group.')
    parser.add_argument('input_file', type=str, help='Input CSV file')
    parser.add_argument('group_file', type=str, help='Group CSV file')
    args = parser.parse_args()
    main(args.input_file, args.group_file)
