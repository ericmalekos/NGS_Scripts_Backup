#!/usr/bin/env bash
###
# Script for single end salmon sequencing. Ignores paired end
# files by excluding _[12].fastq.gz and _R[12].fastq.gz files. Will not ignore 
# PEs that are otherwise formatted
# 
# EX:
#
###

usage() { echo "Usage: $0 [-t <tpm_cutoff>] [-o <output_file>]" 1>&2; exit 0; }

while getopts ":t:o:" arg; do
    case "${arg}" in
        t)
            t=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done


if [ -z "${t}" ] ; then
    t=1
fi


if [ -z "${o}" ] ; then
    o=uniqTscripts.txt
fi

find . -name "quant.sf" | while read FILE; do

        DIR="$(dirname "${FILE}")"
        echo $FILE
        
        awk -v tpm="$t" '$4 > tpm {print $1}' ${FILE} | cut -d"." -f1 >> tscripts.tmp
        sort tscripts.tmp | uniq > ${o}
        rm tscripts.tmp
done
