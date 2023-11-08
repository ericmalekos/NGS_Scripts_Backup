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

echo Running with the following settings:
echo threads \(-t\)          = "${t}"
echo output directory \(-o\) = "${o}"
echo stringency \(-s\)       = "${s}"
echo 

sleep 15s 

find . -not -name "*_[12].fastq.gz" -a -not -name "*_R[12].fastq.gz" -a -name "*.fastq.gz" | while read FILE; do

        DIR="$(dirname "${FILE}")" ; FASTQ="$(basename "${FILE}" .fastq.gz )"
        OUT="${DIR}"/"${o}"
        mkdir -p "${OUT}"
        
        echo Trimming "${FASTQ}"
        
        trim_galore \
        --trim-n  \
        --cores ${t} \
        --stringency ${s} \
        --output_dir "${DIR}"/"${o}"/ \
        "${DIR}"/"${FASTQ}".fastq.gz
        
done
