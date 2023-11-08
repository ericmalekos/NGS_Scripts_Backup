#!/usr/bin/env bash
###
# Script for adapter trimming of PAIRED end reads. Find all _1.fastq.gz in all
# subdirectories and removes adapters based on adapter input. Ignores paired end
# files without "_[12].fastq.gz" suffix and _R[12].fastq.gz files. Will not detect 
# PEs that are otherwise formatted, i.e. use _1.fq.gz, _2.fq.gz NOT _R1.fq.gz, 
# EX:
# ./PE_fastpTrim.sh -a adapter.fa -o FastpQC -t 4
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

find . -name "*_1.fastq.gz" | while read FILE; do

        DIR="$(dirname "${FILE}")" ; FASTQ="$(basename "${FILE}" _1.fastq.gz )"
        OUT="${DIR}"/"${o}"
        mkdir -p "${OUT}"
        
        echo Trimming "${FASTQ}"
        
        fastp \
        --in1 "${DIR}"/"${FASTQ}"_1.fastq.gz \
        --in2 "${DIR}"/"${FASTQ}"_2.fastq.gz \
        --adapter_fasta "${a}" \
        --compression 1 \
        --dont_eval_duplication \
        --thread "${t}" \
        --length_required 19 \
        --json "${OUT}"/"${FASTQ}".json \
        --html "${OUT}"/"${FASTQ}".html \

done
