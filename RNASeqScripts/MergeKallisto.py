#!/usr/bin/env python3.10

import os
import pandas as pd
import argparse

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Merge Kallisto abundance.tsv files and add gene names and types')
parser.add_argument('root_dir', type=str, help='Root directory where all the directories with abundance.tsv files are present')
parser.add_argument('gene_info', type=str, help='TSV file with transcript ID, gene names and gene types')

# Parse command line arguments
args = parser.parse_args()

# Get the root directory and gene_info file from the command-line arguments
root_dir = args.root_dir
gene_info_file = args.gene_info

# Read gene_info_file into a pandas dataframe
gene_info_df = pd.read_csv(gene_info_file, sep="\t", names=["transcript_id", "gene_name", "gene_type"])

# Process transcript_id to remove the .version part
gene_info_df['transcript_id'] = gene_info_df['transcript_id'].str.split(".").str[0]

# Get all subdirectories in the root directory
subdirs = sorted([d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))])

# Loop over all subdirectories
for subdir in subdirs:
    # Create full file path
    file_path = os.path.join(root_dir, subdir, "abundance.tsv")
    
    # Read the tsv file into a pandas dataframe
    df = pd.read_csv(file_path, sep="\t", usecols=["target_id", "tpm"])
    
    # Process target_id to remove the .version part
    df['target_id'] = df['target_id'].str.split(".").str[0]
    
    # Rename the 'tpm' column to 'subdir_tpm' and 'target_id' to 'transcript_id'
    df.rename(columns={"tpm": f"{subdir}_tpm", "target_id": "transcript_id"}, inplace=True)
    
    # Merge df with gene_info_df on transcript_id
    gene_info_df = pd.merge(gene_info_df, df, on="transcript_id", how="left")

# Sort the merged dataframe by 'gene_name'
gene_info_df.sort_values(by='gene_name', inplace=True)

# Print the sorted merged dataframe to the command line
print(gene_info_df.to_csv(sep="\t", index=False))
