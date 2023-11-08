#!/usr/bin/env bash
###
# Script for single end salmon sequencing. Ignores paired end
# files by excluding _[12].fastq.gz and _R[12].fastq.gz files. Will not ignore 
# PEs that are otherwise formatted
# 
# EX:
#  SERNASeqSalmon.sh 
###

usage() { echo "Usage: $0 [-n <index>] [-t <threads>] [-s <read suffix/file extension>] [-o <output_prefix>] [-p <read_prefix>]" 1>&2; exit 0; }

while getopts ":n:o:t:s:r" arg; do
    case "${arg}" in
        n)
            n=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done


if [ -z "${n}" ] || [ -z "${t}" ] ; then
    usage
fi

if [ -z "${s}" ] ; then
    s=.fastq.gz
fi

if [ -z "${p}" ] ; then
    p=trim
fi


if [ -z "${o}" ] ; then
    o=""
    echo No output identifier, saving as 'SAMPLENAME_SalmonQuant'
fi

find . -not -name "*_[12].fastq.gz" -a -not -name "*_R[12].fastq.gz" -a -name "*.fastq.gz" | while read FILE; do

        DIR="$(dirname "${FILE}")" ; SAMPLE="$(basename "${FILE}" "${s}")"

        echo
        echo Sample:  "${SAMPLE1}"

        INDEX="$(basename "${n}")"

        OUT="${DIR}"/"${INDEX}"_"${SAMPLE}"_SalmonQuant/
	    
	    
        salmon quant \
        --index "${n}" \
        --libType A \
        --threads "${t}" \
        --output "${OUT}" \
        -r "${FILE}"
done
