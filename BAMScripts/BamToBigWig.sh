#!/usr/bin/env bash
###
# Finicky and doesnt work well
#
# BamToBigWig -s _merged.bam -t 2
#
###


usage() { echo "Usage: $0 [-p <prefix>] [-s <suffix>] [-t <threads>]" 1>&2; exit 0; }

while getopts ":p:s:t:" arg; do
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
        *)
            usage
            ;;
    esac
done

if [ -z "${s}" ] ; then
    usage
fi

if [ -z "${t}" ] ; then
    t=8
fi

eval "$(conda shell.bash hook)"
conda activate deeptools

find . -name "${p}*${s}" | while read FILE; do

        DIR="$(dirname "$(realpath "${FILE}")")"
        #echo "${DIR}"
        
        DIRBASE="$(basename "${DIR}")"
        
        #echo "${DIRBASE}"
        
        
        SAMPLE="$(basename "${FILE}")"
        SAMPLE=${SAMPLE#"${p}"}
        SAMPLE=${SAMPLE%"${s}"}
        
        echo creating "${SAMPLE}".bw
        
        sleep 10
        
        bamCoverage \
        --bam "${FILE}" \
        --outFileName "${SAMPLE}".bw \
        --numberOfProcessors "${t}" \
        --normalizeUsing CPM \
        --binSize 1
        
done
