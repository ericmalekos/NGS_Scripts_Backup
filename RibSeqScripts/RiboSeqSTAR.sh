#!/usr/bin/env bash
#############
#
# Script for mapping of SINGLE end reads. 
# if z flag is set the input is zipped
#
# EX:
#   RiboSeqSTAR -a ../../../../STAR_Genomes/GencodevM31_index_oh35/ -z -t 16 -s _norRNA.fq.gz
#
#############

# TODO check if its zipped by examining suffix

usage() { echo "Usage: $0 [-a <star_index>] [-z <zipped_input>] [-t <threads>] [-p <prefix>] [-s <suffix>]" 1>&2; exit 0; }

while getopts ":a:z:t:p:s:" arg; do
    case "${arg}" in
        a)
            a=${OPTARG}
            ;;

        t)
            t=${OPTARG}
            ;;
        z)
            z=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${a}" ] || [ -z "${t}" ]; then
    usage
fi

if [ -z "${z}" ] ; then
    z=-
else
    z=zcat
fi


if [ -z "${p}" ] ; then
    p=""
fi

if [ -z "${s}" ] ; then
    s=fastq
fi

find . -name "${p}*${s}" | while read FILE; do

        DIR="$(dirname "${FILE}")"
        
        SAMPLE="$(basename "${FILE}")"
        SAMPLE=${SAMPLE#"$p"}
        SAMPLE=${SAMPLE%"$s"}
        
        OUT="${DIR}"/"${SAMPLE}"_MappedReads
        
        mkdir -p "${OUT}"
        
        echo Trimming and merging "${SAMPLE}"
        echo Output to "${OUT}"
        
        STAR \
        --runMode alignReads \
        --genomeLoad LoadAndKeep \
        --readFilesCommand "${z}" \
        --runThreadN "${t}" \
        --genomeDir  "${a}" \
        --readFilesIn "${FILE}" \
        --outFilterMultimapNmax 1 \
        --outFileNamePrefix "${OUT}"/"${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts TranscriptomeSAM \
        --twopassMode None \
        --outSAMattributes All \
        --outFilterMismatchNmax 3 \
        --alignEndsType EndToEnd \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outMultimapperOrder Random \
        --outSAMmultNmax 1 \
        --limitBAMsortRAM 25000000000

done

STAR \
--runMode alignReads \
--genomeLoad Remove \
--genomeDir  "${a}"
