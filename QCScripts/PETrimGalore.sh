#!/usr/bin/env bash
###
# 
#
#
###

usage() { echo "Usage: $0 [-t <threads>] [-o <output_dir>] [-s <stringency>]" 1>&2; exit 0; }

while getopts ":t:o:s:" arg; do
    case "${arg}" in
        t)
            t=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

#If any input is empty, exit
if [ -z "${t}" ]; then
    t=4
fi

if [ -z "${o}" ]; then
    o=TrimmedReads/
fi

if [ -z "${s}" ]; then
    s=4
fi

find . -name "*_1.fastq.gz" | while read FILE; do

        DIR="$(dirname "${FILE}")" ; FASTQ="$(basename "${FILE}" _1.fastq.gz )"
        
        echo Trimming "${FASTQ}"
        
        trim_galore \
        --trim-n  \
        --cores ${t} \
        --stringency ${s} \
        --output_dir "${DIR}"/"${o}"/ \
        --paired \
        "${DIR}"/"${FASTQ}"_1.fastq.gz \
        "${DIR}"/"${FASTQ}"_2.fastq.gz
        
        mv "${DIR}"/"${o}"/"${FASTQ}"_1_val_1.fq.gz "${DIR}"/"${o}"/trimmed_"${FASTQ}"_1.fastq.gz
        mv "${DIR}"/"${o}"/"${FASTQ}"_2_val_2.fq.gz "${DIR}"/"${o}"/trimmed_"${FASTQ}"_2.fastq.gz
        
done
