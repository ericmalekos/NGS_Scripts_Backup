#!/usr/bin/env bash
###
# USE WITH CAUTION
# This script worked appropriately with the following arguments in the RiboSeq workflow,
# but has not been extensively tested. By piping output to "Merged.out" you can examine
# which BAMs were actually merged.
# MergeBAMReplicates -p "chr_" -s "_Aligned.sortedByCoord.out.bam" -t 2 > Merged.out
#
###
shopt -s expand_aliases
alias samtools='samtools-1.15'
samtools --version | head -1

usage() { echo "Usage: $0 [-p <prefix>] [-s <suffix>] [-t <threads>] [-r <run (y)>]" 1>&2; exit 0; }

while getopts ":p:s:t:r:" arg; do
    case "${arg}" in
        p)
            p=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${p}" ] || [ -z "${s}" ] ; then
    usage
fi

if [ -z "${t}" ] ; then
    t=8
fi

if [ -z "${r}" ] ; then
    find . -name "${p}*[Rr]ep1${s}" | while read FILE; do

            DIR="$(dirname "${FILE}")"
            NOREPDIR=${DIR%[Rr]ep1_MappedReads}
            
            SAMPLE="$(basename "${FILE}")"
            SAMPLE=${SAMPLE#"$p"}
            SAMPLE=${SAMPLE%[Rr]ep1"$s"}
            
            echo merging "${NOREPDIR}"*MappedReads/"$p""${SAMPLE}"[Rr]ep[123456789]"${s}"            

    done
    echo
    echo If this looks correct, run again with \`-r\` flag set
    exit
fi        

find . -name "${p}*[Rr]ep1${s}" | while read FILE; do

        DIR="$(dirname "${FILE}")"
        NOREPDIR=${DIR%[Rr]ep1_MappedReads}
        
        SAMPLE="$(basename "${FILE}")"
        SAMPLE=${SAMPLE#"$p"}
        SAMPLE=${SAMPLE%[Rr]ep1"$s"}
        
        echo merging "${NOREPDIR}"*MappedReads/"$p""${SAMPLE}"[Rr]ep[123456789]"${s}"            

        samtools merge -f -@"${t}" "${NOREPDIR}"merged.bam "${NOREPDIR}"*MappedReads/"$p""${SAMPLE}"[Rr]ep[123456789]"${s}"
        samtools index -@"${t}" "${NOREPDIR}"merged.bam

done
