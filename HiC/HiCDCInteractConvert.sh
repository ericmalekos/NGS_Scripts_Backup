#!/usr/bin/env bash

####################################################################
####################################################################
#
# Convert HiCDC+ output to bedInteract for viewing in genome browser
# and an interact file for viewing with pygenometracks
#
####################################################################
####################################################################

# Check if the file argument is passed
if [ $# -eq 0 ]; then
    echo "Please provide the file name as an argument."
    exit 1
fi

# Store the file name in a variable
FILE=$1
NAME=${FILE%.*}
# Check if the file exists
if [ ! -f ${FILE} ]; then
    echo "File ${FILE} not found."
    exit 1
fi

# Rearrange the columns using awk
# Expected input format
# chrI	startI	endI	chrJ	startJ	endJ	Distance	counts	pvalue	qvalue	mu	sdev
# 
# Desired
# chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand

# $7 != 0 to ignore self-self interactions
tail -n +2 ${FILE} | awk -v OFS='\t' ' $7 != 0 { print $1, $2, $6, " ", $8, "0", ".", "0", $1, $2, $3, " ", ".", $4, $5, $6, " ", "."}' > ${NAME}.bedInteract
tail -n +2 ${FILE} | awk -v OFS='\t' ' $7 != 0 { print $1, $2, $3, $4, $5, $6, $8}' > ${NAME}.pairs

# Confirm the rearrangement
echo "Complete"


