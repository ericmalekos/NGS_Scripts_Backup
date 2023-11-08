#!/usr/bin/env bash

####################################################################
####################################################################
#
# Removes everything that is non-canonical from each bam. 
# Optionally removes overlapping entries in a given bed/gtf file
#
# Made for use with STAR output in mind "Aligned.sortedByCoord.out.bam"
# Works with either human or mouse, but results in annoying warning for mouse data:
# "[main_samview] region "chr20" specifies an invalid region or unknown reference. Continue anyway."
#
# With -d flag set to "y", the original bam and bam.bai/csi will be deleted 
# following samtools quickcheck on the new bam
# USE `-d y` WITH CARE
#
# Use -i to remove reads from unwanted positions, -i snoRNAs.bed
# same strand and 33% overlap for removal hardcoded. 
####################################################################
####################################################################
shopt -s expand_aliases
alias samtools='samtools-1.15'
samtools --version | head -1

usage() { echo "Usage: $0 [-t <threads>] [-d <delete_old_bams>] [-i <remove_overlap_reads_bed/gtf>]" 1>&2; exit 0; }

while getopts ":t:d:i:" arg; do
    case "${arg}" in
        t)
            t=${OPTARG}
            ;;
        d)
 	        d=${OPTARG}
            ;;
        i)
	        i=${OPTARG}
            ;;
        
        *)
            usage
            ;;
    esac
done

if [ "${t}" == "" ]; then
  t=1
fi	

####################################################################


find . -name "*Aligned.sortedByCoord.out.bam" | while read FILE; do
    
    echo "${FILE}"
    
    DIR="$(dirname "${FILE}")"
    BAM="$(basename "${FILE}")"
    
    #index bam if it's not index
    
    [ ! -f "${FILE}".bai -a ! -f "${FILE}".csi ] && \
     samtools index -b -@"${t}" "${FILE}"
    
    # If bed/gtf/etc. passed, then remove intersecting reads
    # bedtools intersect on same strand, minimum 33% of read overlaps unwanted bed entry
    # -v flags 
    if [ -z "${i}" ] ; then
    	samtools view -bh -@"${t}" "${FILE}" chr{1..22} chrX chrY chrM |\
    	samtools sort -l1 -@"${t}" - > "${DIR}"/chr_"${BAM}"
    	if [ $? -eq 1 ] ; then
    	  exit 1
    	fi
    else
    	samtools view -bh -@"${t}" "${FILE}" chr{1..22} chrX chrY chrM | \
    	samtools sort -l1 -@"${t}" - | \
    	bedtools intersect  -v -s -f 0.33 -a stdin -b "${i}" > "${DIR}"/chr_"${BAM}"
    	if [ $? -eq 1 ] ; then
    	  exit 1
    	fi
    fi
    
    samtools index -@"${t}" "${DIR}"/chr_"${BAM}"
    
    # If a new bam file was generated and -d is set, delete the original
    samtools quickcheck "${DIR}"/chr_"${BAM}"
    if [ $? -eq 0 ] && [ "${d}" == "y" ]; then
    	rm "${FILE}"
    	rm "${FILE}""."[bc][as]"i"
    fi
    
done

