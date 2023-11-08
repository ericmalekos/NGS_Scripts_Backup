#!/usr/bin/env bash
###
#
#   ParallelRiboTricer.sh -i  ./RiboTricerIndex/AgatMergeGen31LNCRNASmale_candidate_orfs.tsv 
#   -o ./RiboTricerPredictions/ -t 10 -p 0.418 -s _merged.bam
#
# 
###

usage() { echo "Usage: $0 [-i <input_directory>] [-o <output_directory>] [-t <threads>] [-s <suffix>] [-p phase_cutoff_score]" 1>&2; exit 0; }

while getopts ":i:o:t:s:p:" arg; do
    case "${arg}" in
        i)
            i=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

#if [ -z "${i}" ] || [ -z "${o}" ] || [ -z "${t}" ] || [ -z "${p}" ] || [ -z "${s}" ] ; then
#    usage
#fi

echo "Index:  ${i}"
echo "Threads:        ${t}"
echo "Output Directory: ${o}"
echo "phasecutoff:        ${p}"
echo "Suffix: ${s}"

# semaphore code from https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop/103922
# initialize a semaphore with a given number of tokens - set with the `-t` parameter
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}


RiboTricer () {
    local BAM="${1}"
    local DIR="$(dirname "${BAM}")" ;
    local PREFIX="$(basename "${BAM}" "${s}")" ; 
    mkdir -p ${o}/"${DIR}"_"${PREFIX}"/ ; 
    ribotricer detect-orfs \
      --bam ${BAM} \
      --ribotricer_index ${i} \
      --prefix "${DIR}"_"${PREFIX}"/"${PREFIX}" \
      --read_lengths 27,28,29,30,31,32,33 \
      --phase_score_cutoff ${p} ; 
}


open_sem $t
find . -name *$s -a -not -name "*LTM*"  -a -not -name "*HARR*"| while read FILE; do \
    echo ${FILE}
    run_with_lock RiboTricer ${FILE}
done

