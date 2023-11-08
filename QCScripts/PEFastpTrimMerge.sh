#!/usr/bin/env bash
###
# Script for adapter trimming and MERGING of PAIRED end reads. Find all _1.fastq.gz in all
# subdirectories and removes adapters based on adapter input. Ignores paired end
# files without "_[12].fastq.gz" suffix and _R[12].fastq.gz files. Will not detect 
# PEs that are otherwise formatted, i.e. use _1.fq.gz, _2.fq.gz NOT _R1.fq.gz, 
# EX:
# ./PE_fastpTrimMerge.sh -a adapter.fa -o FastpQC -t 4 -l 15 -s fastq.gz
###

usage() { echo "Usage: $0 [-a <adapter_fasta>] [-o <output_directory>] [-t <threads>] [-l <overlap_length>] [-s suffix]" 1>&2; exit 0; }

while getopts ":a:o:t:l:s:" arg; do
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
        l)
            l=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

#If any input is empty, exit
if [ -z "${a}" ] || [ -z "${o}" ] || [ -z "${t}" ] || [ -z "${l}" ] || [ -z "${s}" ] ; then
    usage
fi

find . -name "*_1""${s}" | while read FILE; do

        DIR="$(dirname "${FILE}")" ; FASTQ="$(basename "${FILE}" _1"${s}")"
        OUT="${DIR}"/"${o}"
        mkdir -p "${OUT}"
        
        echo Trimming and merging "${FASTQ}"
        
        fastp \
        --in1 "${DIR}"/"${FASTQ}"_1"${s}" \
        --in2 "${DIR}"/"${FASTQ}"_2"${s}" \
        --adapter_fasta "${a}" \
        --compression 1 \
        --dont_eval_duplication \
        --thread "${t}" \
        --length_required 19 \
        --json "${OUT}"/"${FASTQ}".json \
        --html "${OUT}"/"${FASTQ}".html \
        --merge \
        --merged_out "${DIR}"/trimmed_merged_"${FASTQ}".fastq.gz \
        --overlap_len_require "${l}"

done
