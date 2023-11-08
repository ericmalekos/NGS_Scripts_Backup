#!/usr/bin/env bash
###
# Script for adapter trimming of SINGLE end reads. Find all fastq.gz in all
# subdirectories and removes adapters based on adapter input. Ignores paired end
# files by excluding _[12].fastq.gz and _R[12].fastq.gz files. Will not ignore 
# PEs that are otherwise formatted
# 
# EX:
#  SEFastpTrim -a adapter.fa -o FastpQC -t 4 
###

usage() { echo "Usage: $0 [-a <adapter_fasta>] [-o <output_directory>] [-t <threads>]" 1>&2; exit 0; }

while getopts ":a:o:t:" arg; do
    case "${arg}" in
        a)
            a=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

#If any input is empty, exit
if [ -z "${a}" ] || [ -z "${o}" ] || [ -z "${t}" ]; then
    usage
fi

find . -not -name "*_[12].fastq.gz" -a -not -name "*_R[12].fastq.gz" -a -name "*.fastq.gz" | while read FILE; do

        DIR="$(dirname "${FILE}")" ; FASTQ="$(basename "${FILE}" .fastq.gz )"
        OUT="${DIR}"/"${o}"
        mkdir -p "${OUT}"
        
        echo Trimming and merging "${FASTQ}"
        
        fastp \
        --in1 "${FILE}" \
        --out1 "${DIR}"/trimmed_"${FASTQ}".fastq.gz \
        --adapter_fasta "${a}" \
        --compression 1 \
        --dont_eval_duplication \
        --thread "${t}" \
        --length_required 19 \
        --json "${OUT}"/"${FASTQ}".json \
        --html "${OUT}"/"${FASTQ}".html
done
