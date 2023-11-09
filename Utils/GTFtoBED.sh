#!/usr/bin/env bash

####################################################################
####################################################################
# 
# convert gtf/gff to bed. Checks the file extension to determine format
#
# pass input file as first argument, file of transcript IDs as optional second argument
# 
# Only outputs lines with "transcript" or "RNA" in third column (i.e. no exons)
# Output: <chrom> <pos1> <pos2> <tscriptID> <1> <+/->
#
####################################################################
####################################################################



if [ -z "$1" ] ; then
	>&2 echo
	>&2 echo "Provide input GTF/GFF"
	>&2 echo "Optionally provide a line separated list of transcript IDs to subset"
	>&2 echo "Output prints to stdout"
	>&2 echo
	exit
elif [ "${1: -4}" == ".gtf"  ] && [ "$2" ] ; then
	>&2 echo
	>&2 echo "GTF detected, piping through gffread"
	>&2 echo "Transcript ID file detected, subsetting by ID"
	>&2 echo
	gffread -E "$1" -o- | grep -F -f "$2" - | \
	awk '$3 == "RNA" || $3 == "transcript" {print}' | cut -f1 -d";" | awk 'BEGIN {OFS="\t"};{gsub (/ID=/,"",$9) ; print $1,$4,$5,$9,1,$7}' \
	| bedtools sort -i stdin

elif  [ "${1: -4}" == ".gtf" ] && [ -z "$2" ] ; then
	>&2 echo
	>&2 echo "GTF detected, piping through gffread"
	>&2 echo "No Transcript ID file detected, outputing all results"
	>&2 echo
	gffread -E "$1" -o- | \
	awk '$3 == "RNA" || $3 == "transcript" {print}' | cut -f1 -d";" | awk 'BEGIN {OFS="\t"};{gsub (/ID=/,"",$9) ; print $1,$4,$5,$9,1,$7}' \
	| bedtools sort -i stdin
	
elif [ "${1: -4}" == ".gff" ] && [ "$2" ] ; then
	>&2 echo
	>&2 echo "Transcript ID file detected, subsetting by ID"
	>&2 echo
	grep -F -f "$2" "$1" | \
	awk '$3 == "RNA" || $3 == "transcript" {print}' | cut -f1 -d";" | awk 'BEGIN {OFS="\t"};{gsub (/ID=/,"",$9) ; print $1,$4,$5,$9,1,$7}' \
	| bedtools sort -i stdin


elif  [ "${1: -4}" == ".gff"  ] && [ -z "$2" ] ; then
	>&2 echo
	>&2 echo "No Transcript ID file detected, outputing all results"
	>&2 echo
	cat "$1" | \
	awk '$3 == "RNA" || $3 == "transcript" {print}' | cut -f1 -d";" | awk 'BEGIN {OFS="\t"};{gsub (/ID=/,"",$9) ; print $1,$4,$5,$9,1,$7}' \
	| bedtools sort -i stdin
else 
	echo
	echo "Check input files"
	echo "Expected input:"
	echo "GTFtoBED.sh <input.gtf/gff> <tscriptIDs.txt>"
	echo
fi
