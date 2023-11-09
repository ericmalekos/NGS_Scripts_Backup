#!/bin/bash

# Check if input and gff file arguments were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <gff_file>"
    exit 1
fi

# Get the input file from command-line argument
input_file="$1"
gff_file="$2"

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
}' "$gff_file" | sort -k1,1 > gff_mapping.tsv

# Now, process the input file
awk '
BEGIN {FS=OFS="\t"}
{
    pair = $3 OFS $4;
    count[pair]++;
}
END {
    for (pair in count) {
        split(pair, a, OFS);
        split(a[2], b, "_");
        transcript = b[1];
        for (i=2; i<=length(b)-3; i++)
            transcript = transcript "_" b[i];
        split(transcript, t, ".");
        print a[1], a[2], count[pair], t[1]
    }
}' "$input_file" | sort -t $'\t' -k4,4 > result.tsv

# Join the result with the gff mapping
join -t $'\t' -1 4 -2 1 result.tsv gff_mapping.tsv | sort -t $'\t' -k6,6 -k5,5 -k4,4r | cut -f2- | cat

# Cleanup intermediate files
#rm gff_mapping.tsv result.tsv
