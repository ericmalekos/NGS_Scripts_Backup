#!/usr/bin/env bash

####################################################################
####################################################################
# 
# determine likely crispr cut site by shifting by 3 nts
#
# pass input file as argument
# 
# Output shifted bed -3 on the positive strand, +3 on the negative strand
#
####################################################################
####################################################################



if [ -z "$1" ] || [ "${1: -4}" == ".bed"  ] ; then
	>&2 echo
	>&2 echo "Provide input BED"
	>&2 echo "Output prints to stdout"
	>&2 echo
	exit
fi


awk -v OFS='\t' 'BEGIN { if ( $6 == "-" ) print}' "$1"

