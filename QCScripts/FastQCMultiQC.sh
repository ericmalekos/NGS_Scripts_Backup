#!/usr/bin/env bash

####################################################################
####################################################################
# This script runs fastqc and multiqc on all fastq files
# directory structure is assumed to look like:
#
# in "/RunFromDir" execute with
#   `FastQCMultiQC -t 12 -o RawQC -s fq.gz`
####################################################################
####################################################################

usage() { echo "Usage: $0 [-t <threads>] [-o <output_directory_name>] [-s <suffix>] -d [single_directory]" 1>&2; exit 0; }

while getopts ":t:o:s:d:" arg; do
    case "${arg}" in
        t)
            t=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
	    s)
            s=${OPTARG}
            ;;
	    d)
            d=${OPTARG}
            ;;

        *)
            usage
            ;;
    esac
done

if [ -z ${t} ] || [ -z ${o} ] ; then
	usage
fi

if [ -z ${s} ] ; then
	s=fastq.gz
	echo Assuming fastq.gz file extension \( change with \`-s\` parameter \)
fi

if [ -z ${d} ] ; then
	d=.
	echo Searching all subdirectories for files ending in "${s}"
else
    echo Searching in "${d}" and its subdirectories for files ending in "${s}"
fi


####################################################################
# QC on raw reads
# for each subdirectory with fastq.gz files, perform fastqc and multiqc
# output:
#    Directory with fastqc and multiqc reports

eval "$(conda shell.bash hook)"
conda activate multiqc

for DIR in $(find ${d} -name "*${s}" -printf '%h\n' | sort -u); do
    echo "${DIR}"
    echo Output to "${DIR}"/"${o}"
    mkdir -p "${DIR}"/"${o}"
    fastqc --threads "${t}" --outdir "${DIR}"/"${o}"  "${DIR}"/*"${s}"
    multiqc --outdir "${DIR}"/"${o}" .
done

conda deactivate
