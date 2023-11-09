#!/bin/bash

# Check if input and gff file arguments were provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <gff_file>"
    exit 1
fi

# Get the input file from command-line argument
gff_file="$1"

# First, process the GFF file to extract the transcript_id, gene_name and gene_type
awk '
BEGIN {FS=OFS="\t"}
$3 == "transcript" {
    split($9, a, ";");
    for (i in a) {
        split(a[i], b, "=");
        if (b[1] == "transcript_id") {
            split(b[2], t, ".");
            transcript_id = t[1];
        }
        else if (b[1] == "gene_name") {
            gene_name = b[2];
        }
        else if (b[1] == "gene_type") {
            gene_type = b[2];
        }
    }
    print transcript_id, gene_name, gene_type
}
$3 == "RNA" {
    split($9, a, ";");
    for (i in a) {
        split(a[i], b, "=");
        if (b[1] == "transcript_id") {
            split(b[2], t, ".");
            transcript_id = t[1];
        }
        else if (b[1] == "gene_id") {
            gene_name = b[2];
        }
    }
    gene_type = $1 ":" $4 "-" $5;
    print transcript_id, gene_name, gene_type
}' "$gff_file" | sort -k1,1 | uniq | cat
