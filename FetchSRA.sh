#!/usr/bin/env bash
###
# Script for downloading a list of SRA files.
# 
# EX:
#  FetchSRA -o DATA/ -s SRA/SRARunTable.txt
###

usage() { echo "Usage: $0 [-s <sra_table>] [-o <output_directory>] [-t <threads>]" 1>&2; exit 0; }

while getopts ":s:o:t:" arg; do
    case "${arg}" in
        s)
            s=${OPTARG}
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

if [ -z "${s}" ] || [ -z "${o}" ]; then
    usage
fi

if [ -z "${t}" ] ; then
    t=4
fi

sed 1d ${s} | while IFS=, read -r SRA unwanted
do
	echo "Retrieving $SRA"

	/public/groups/carpenterlab/people/emalekos/bin/sratoolkit.2.11.0/bin/fasterq-dump.2.11.0 \
		$SRA \
		--outdir ${o} \
		--threads ${t} \
		--details \
		--progress
	# zip new fastqs
	 
	for FASTQ in ${o}/*fastq; do pigz --fast --processes ${t} ${FASTQ} ; done 
done
