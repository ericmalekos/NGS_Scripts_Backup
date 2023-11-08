#!/usr/bin/env bash
###
# Created: 08-03-2021
# Updated: 09-12-2021
#
# Script for Mapping STAR
# EX:
# PERNASeqSTAR \
#   -n /public/groups/carpenterlab/people/emalekos/STAR_Genomes/merged_nanopore_gencode25_index_oh75 \
#   -i /public/groups/carpenterlab/people/emalekos/smORFs/DATA/Flavell/raw_reads/ \
#   -t 4 \
#   -o /public/groups/carpenterlab/people/emalekos/smORFs/DATA/Flavell/
#   -s fq.gz
#
###

usage() { echo "Usage: $0 [-i <input_directory>] [-n <index>] [-t <threads>] [-s <suffix/file extension>]" 1>&2; exit 0; }

while getopts ":i:n:o:t:s:" arg; do
    case "${arg}" in
        i)
            i=${OPTARG}
            ;;
        n)
            n=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done


if [ -z "${i}" ] || [-z "${n}" ] || [-z "${t}" ] || [ -z "${s}" ]; then
    usage
fi

echo "Input Directory:  ${i}"
echo "STAR Index:       ${n}"
echo "Output Directory: ${o}"
echo "Threads:          ${t}"

#for paired end reads

ext=fastq.gz

for READ in ${i}/*_1."${s}"
do
        SAMPLE="$(basename "${READ}" _1."${s}")"
        SAMPLE1="${SAMPLE}"_1."${s}"
        SAMPLE2="${SAMPLE}"_2."${s}"
        echo
        echo First Read:  "${SAMPLE1}"
        echo Second Read: "${SAMPLE2}"

        OUT="${i}"/"${SAMPLE}"_MappedReads
	    
	    mkdir -p "${OUT}"
	    

        STAR \
        --genomeDir ${n} \
        --readFilesIn \
        ${i}/"${SAMPLE1}" \
        ${i}/"${SAMPLE2}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All \
        --twopassMode Basic \
        --outFilterMultimapNmax 10 \
        --quantMode TranscriptomeSAM GeneCounts \
        --runThreadN ${t} \
        --alignEndsType Local \
        --outFileNamePrefix "${OUT}"/${SAMPLE}  
done
