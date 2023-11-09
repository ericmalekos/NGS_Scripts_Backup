#!/usr/bin/env bash
###
#
#
# Ribotish quality command on CHX reads
#
# EX:
#   RibotishQC -g ../../../../annotations/mouse/M31/CCDStagged_gencode.vM31.annotation.gtf -t 24
#   
# NOTE Ribotish does not exit easily with CTRL+c. Instead stop job with CTRL+z
# get job id with   `jobs -p` 
# end job with      `kill -9 $(jobs -p)`
###

usage() { echo "Usage: $0 [-g <gtf_with_CDS>] [-t <threads>] [-s <shortest_read_length>] [-l <longest_read_length>] [-u suffix]" 1>&2; exit 0; }

while getopts ":g:t:s:l:u:" arg; do
    case "${arg}" in
        g)
            g=${OPTARG}
            ;;

        t)
            t=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        u)
            u=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# If any input is empty, exit
if [ -z "${g}" ] || [ -z "${t}" ]; then
    usage
fi

if [ -z "${s}" ] ; then
    s=19
fi

if [ -z "${l}" ] ; then
    l=37
fi

if [ -z "${u}" ] ; then
    u="Aligned.sortedByCoord.out.bam"
fi

eval "$(conda shell.bash hook)"
conda activate ribotish_python2_7


find . -name "*""${u}" | while read FILE; do

        DIR="$(dirname "${FILE}")"
        
        #SAMPLE="$(basename "${FILE}")"
        #SAMPLE=${SAMPLE#"$p"}
        #SAMPLE=${SAMPLE%"$s"}
        
        OUT="${DIR}"/RibotishQC/
        
        mkdir -p "${OUT}"
        
        echo Output to "${OUT}"
        
        ribotish quality \
        -b "${FILE}" \
        -g "${g}" \
        -f "${OUT}"/qual.pdf \
        -r "${OUT}"/para.py \
        -p "${t}" \
        --tis \
        -l "${s}","${l}"

done
