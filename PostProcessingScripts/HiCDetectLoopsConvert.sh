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
#chr1	x1	x2	chr2	y1	y2	pvalue
# 
# Desired for bigInteract
# chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand

awk -v OFS='\t' '{ print $1, $2, $6, ".", 1, "0", ".", "0", $1, $2, $3, " ", ".", $4, $5, $6, " ", "."}' ${FILE} > ${NAME}.bedInteract

awk -v OFS='\t' ' $7 != 0 { print $1, $2, $3, $4, $5, $6, 1 }' ${FILE} > ${NAME}.pairs \


# Confirm the rearrangement
echo "Complete"


