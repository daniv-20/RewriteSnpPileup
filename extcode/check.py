import gzip
import csv
import os
import pandas as pd
import getplpcsv as gc
import argparse

# Print the current working directory
print("Current working directory:", os.getcwd())

## write script to compare the two results

def modify_filename(file_name):
    # Remove '.gz' from the end if it exists
    if file_name.endswith('.gz'):
        file_name = file_name[:-3]
    
    # Add '.csv' if it doesn't already end with '.csv'
    if not file_name.endswith('.csv'):
        file_name += '.csv'
    
    return file_name

def compare(python_out_gz, r_out_gz):
    python_out = modify_filename(python_out_gz)
    r_out = modify_filename(r_out_gz)

## get csvs
    gc.convert_gzip_to_csv(python_out_gz, python_out)
    gc.convert_gzip_to_csv(python_out_gz, python_out)

    py = pd.read_csv(python_out, low_memory = False)
    r = pd.read_csv(r_out, low_memory = False)

    ## for each position in each chromosome, if the position is in both files, 
    # check that the results are the same
    compare_df = pd.merge(py, r, on=["Chromosome", "Position"], how = 'outer', suffixes=('_py', '_r'))

    # Assuming compare_df has already been created
    compare_df['sus'] = 0  # Initialize 'sus' column to 0

    # Loop through all columns from df1 and check corresponding columns from df2
    for col in py.columns:
        if col not in ['Chromosome', 'Position']:
            col_py = f"{col}_py"
            col_r = f"{col}_r"
            if col_py in compare_df.columns and col_r in compare_df.columns:
                # Check if values are different and update 'sus'
                compare_df['sus'] |= (compare_df[col_py] != compare_df[col_r]).astype(int)
    return compare_df

def results_summary(compare_df):
    agree = compare_df.loc[compare_df.sus == 0, ]
    print("Agree: ", agree.shape[0])
    disagree = compare_df.loc[compare_df.sus == 1, ]
    print("Disagree: ", disagree.shape[0])

def main():
    parser = argparse.ArgumentParser(description="Check Results")
    parser.add_argument("-py", "--python", required=True, help="Path to the python output file.")
    parser.add_argument("-r", "--r", required=True, help="Path to the r output file.")
    args = parser.parse_args()
    compare_df = compare(args.python, args.r)
    results_summary(compare_df)

if __name__ == "__main__":
    main()