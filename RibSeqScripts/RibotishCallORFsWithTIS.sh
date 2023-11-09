#!/usr/bin/env bash
###
# First run to determine if the correct bams are being identified
#   RibotishCallORFs -g genome.fa -f annotation.gtf -t 4 -p "." -s .para.py -b _merged.bam -i RibotishQCDirectory
# If correct, run again with the `d` flag set
#
# NOTE Ribotish does not exit easily with CTRL+c. Instead stop job with CTRL+z
# get job id with   `jobs -p` 
# end job with      `kill -9 $(jobs -p)`
#
###

usage() { echo "Usage: $0 [-g <genome>] [-f <gtf>] [-t <threads>] [-p <prefix_seperator>] [-s <offset_suffix>] [-b <bam_suffix>] [-o <output_dir>] [-i <input_dir>] [-m tis_bam] [-e tis_offset] [-h harr] [-d <proceed_with_ORFCalling>]" 1>&2; exit 0; }

while getopts ":g:f:t:p:s:b:o:i:d:m:e:h:" arg; do
    case "${arg}" in
        g)
            g=${OPTARG}
            ;;
        f)
            f=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        d)
            d=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        e)
            e=${OPTARG}
            ;;
        h)
            h=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${g}" ] || [ -z "${f}" ] || [ -z "${t}" ] || [ -z "${i}" ] || [ -z "${e}" ] || [ -z "${m}" ] ; then
    usage
fi

if [ -z "${p}" ] ; then
    p="."
fi

if [ -z "${o}" ] ; then
    o=./ORFsOut
fi

if [ -z "${s}" ] ; then
    s=.para.py
fi

if [ -z "${b}" ] ; then
    b=_merged.bam
fi


if [ -z "${d}" ] ; then
    d=n
fi

if [ -z "${h}" ] ; then
    h=n
fi

# Check that the QC files have corresponding bam files

find "${i}" -name *"${s}" | while read FILE ; do

        echo Input QC offset file: "${FILE}"
        
        TOPDIR="$(basename "${FILE}")"
        TOPDIR=${TOPDIR%%"${p}"*}  
        
        BAM="$(basename "${FILE}" "${s}")"
        BAM="${BAM#*"${p}"}"
        BAM="${BAM}""${b}"
        
        echo Bam file "${BAM}"
        
        BAMPATH=$(find "${TOPDIR}" -name "${BAM}")
        
        if test -f "${BAMPATH}"; then
            echo "${BAMPATH}" exists
        else
            echo
            echo COULD NOT FIND "${BAMPATH}"
            exit
        fi    
done

# Now start ribotish

if [ "${d}" != "n" ] ; then

    echo
    echo continuing with ORF Calling
    echo


    eval "$(conda shell.bash hook)"
    conda activate ribotish_python2_7

    find "${i}" -name *"${s}" | while read FILE; do

            #echo Input QC offset file: "${FILE}"
            #exit
            TOPDIR="$(basename "${FILE}")"
            TOPDIR=${TOPDIR%%"${p}"*}  
            
            #echo Directory "${TOPDIR}"
            
            BAM="$(basename "${FILE}" "${s}")"
            BAM="${BAM#*"${p}"}"
            BAM="${BAM}""${b}"
            
            BAMPATH=$(find "${TOPDIR}" -name "${BAM}")
            
            OUTDIR="${o}"/"${TOPDIR}"
            mkdir -p "${OUTDIR}"
            
            echo Offset file "${FILE}"
            echo BAM "${BAMPATH}"
            
            if [ "${d}" != "n" ] ; then 
            echo
            echo "\`-h\` not set, TIS input is not harringtonine"
            echo 
                ribotish predict \
                    -t "${m}" \
                    --tispara "${e}" \
                    -b "${BAMPATH}" \
                    -g "${g}" \
                    -f "${f}" \
                    -o "${OUTDIR}"/"${BAM}"ORFs.tsv \
                    --ribopara "${FILE}" \
                    --geneformat gtf \
                    --alt \
                    --longest \
                    --minaalen 10 \
                    --inframecount \
                    --fpth 1 \
                    --allresult "${OUTDIR}"/"${BAM}"AllORFs.tsv \
                    --transprofile "${OUTDIR}"/"${BAM}".transprofile \
                    -p "${t}" \
                    --seq \
                    --aaseq \
                    --blocks
           else
            echo
            echo "\`-h\` set, TIS input is harringtonine"
            echo 
                ribotish predict \
                    --harr \
                    -t "${m}" \
                    --tispara "${e}" \
                    -b "${BAMPATH}" \
                    -g "${g}" \
                    -f "${f}" \
                    -o "${OUTDIR}"/"${BAM}"ORFs.tsv \
                    --ribopara "${FILE}" \
                    --geneformat gtf \
                    --alt \
                    --longest \
                    --minaalen 10 \
                    --inframecount \
                    --fpth 1 \
                    --allresult "${OUTDIR}"/"${BAM}"AllORFs.tsv \
                    --transprofile "${OUTDIR}"/"${BAM}".transprofile \
                    -p "${t}" \
                    --seq \
                    --aaseq \
                    --blocks
          fi
         
    done
fi

