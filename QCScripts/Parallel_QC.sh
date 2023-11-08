#!/usr/bin/env bash
###
# Created: 08-05-2021
# Updated: 08-05-2021
#
# Script for Trimming SE Reads in parallel
#
# TODO: make output directory work
###

usage() { echo "Usage: $0 [-i <input_directory>] [-o <output_directory>] [-p <processes>]" 1>&2; exit 0; }

# semaphore code taken from https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop/103922
# initialize a semaphore with a given number of tokens - set with the `-p` parameter
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




while getopts ":i:o:p:" arg; do
    case "${arg}" in
        i)
            i=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${i}" ] || -z "${n}" ] || [ -z "${o}" ]; then
    usage
fi

echo "Input Directory:  ${i}"
echo "Processes:        ${p}"
echo "Output Directory: ${o}"


if [ -d ${o} ]
then
    echo "Directory ${o} exists."
else
    echo "Directory ${o} does not exists, creating directory."
    mkdir ${o}
fi

# Example
# fqc () {
#     local reads=$1
#     /public/groups/carpenterlab/people/emalekos/bin/FastQC/fastqc \
#     ${i}/${reads}.fastq.gz -o ${o}
# }
#
# open_sem $p
# for reads in ${i}/*fastq.gz; do
#     echo "$(basename "$reads" .fastq.gz)"
#     run_with_lock trim "$(basename "$reads" .fastq.gz)"
# done
