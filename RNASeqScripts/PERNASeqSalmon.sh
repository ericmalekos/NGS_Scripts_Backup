#!/usr/bin/env bash
###
#
# Script for quantifyin with Salmon
# EX:
#   PERNASeqSalmon.sh -n ../../Salmon_Index/Mm39Gencode31/gencode.vM31.comprehensive.index/ -t 8
#
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
    s=fastq.gz
fi

if [ -z "${p}" ] ; then
    p=trim
fi


if [ -z "${o}" ] ; then
    o=""
    echo No output identifier, saving as 'SAMPLENAME_SalmonQuant'
fi

echo "Salmon Index:     ${n}"
echo "Output Directory: ${o}"
echo "Threads:          ${t}"

find . -name "${p}*"_1."${s}" | while read FILE ; do

        echo $FILE

        DIR="$(dirname "${FILE}")" ; SAMPLE="$(basename "${FILE}" _1."${s}")"
        SAMPLE1="${SAMPLE}"_1."${s}"
        SAMPLE2="${SAMPLE}"_2."${s}"
        echo
        echo First Read:  "${SAMPLE1}"
        echo Second Read: "${SAMPLE2}"

        INDEX="$(basename "${n}")"

        OUT="${DIR}"/"${INDEX}"_"${SAMPLE}"_SalmonQuant/
	    
        salmon quant \
        --index "${n}" \
        --libType A \
        --threads "${t}" \
        --output "${OUT}" \
        --mates1 "${DIR}"/"${SAMPLE}"_1."${s}" \
        --mates2 "${DIR}"/"${SAMPLE}"_2."${s}"
        
done
