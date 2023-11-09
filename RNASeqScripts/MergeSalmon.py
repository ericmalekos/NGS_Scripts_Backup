#!/usr/bin/env python3.10
import os
import pandas as pd
import argparse

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Merge Salmon quant.sf files and add gene names and types')
parser.add_argument('root_dir', type=str, help='Root directory where all the directories with quant.sf files are present')
parser.add_argument('gene_info', type=str, help='TSV file with transcript ID, gene names and gene types')

# Parse command line arguments
args = parser.parse_args()

# Get the root directory, gene_info file from the command-line arguments
root_dir = args.root_dir
gene_info_file = args.gene_info

# Read gene_info_file into a pandas dataframe
gene_info_df = pd.read_csv(gene_info_file, sep="\t", names=["transcript_id", "gene_name", "gene_type"])

# Process transcript_id to remove the .version part
gene_info_df['transcript_id'] = gene_info_df['transcript_id'].str.split(".").str[0]

# Collect directories containing quant.sf files
subdirs = [os.path.dirname(os.path.join(root, name)) 
           for root, dirs, files in os.walk(root_dir) 
           for name in files 
           if name == 'quant.sf']

# Sort the subdirs list
subdirs.sort()

# Loop over all sorted subdirectories
for subdir in subdirs:
    # Create full file path
    file_path = os.path.join(subdir, "quant.sf")
    
    # Read the tsv file into a pandas dataframe
    df = pd.read_csv(file_path, sep="\t", usecols=["Name", "TPM"])
    
    # Process Name to remove the .version part
    df['Name'] = df['Name'].str.split(".").str[0]
    
    # Rename the 'TPM' column to 'subdir_TPM' and 'Name' to 'transcript_id'
    df.rename(columns={"TPM": f"{os.path.basename(subdir)}_TPM", "Name": "transcript_id"}, inplace=True)
    
    # Merge df with gene_info_df on transcript_id
    gene_info_df = pd.merge(gene_info_df, df, on="transcript_id", how="left")

# Sort the merged dataframe by 'gene_name'
gene_info_df.sort_values(by='gene_name', inplace=True)

# Specify columns to keep at the front
first_columns = ['transcript_id', 'gene_name', 'gene_type']

# Get the remaining columns and sort them
remaining_columns = gene_info_df.columns.tolist()
remaining_columns = [col for col in remaining_columns if col not in first_columns]

# Combine the first and remaining columns and rearrange dataframe
final_columns_order = first_columns + remaining_columns
gene_info_df = gene_info_df[final_columns_order]

# Print the sorted merged dataframe to the command line
print(gene_info_df.to_csv(sep="\t", index=False))
