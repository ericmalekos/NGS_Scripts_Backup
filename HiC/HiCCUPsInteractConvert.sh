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
#chr1	x1	x2	chr2	y1	y2	name	score	strand1	strand2	color	observed	expectedBL	expectedDonut	expectedH	expectedV	fdrBL	fdrDonut	fdrH	fdrV	numCollapsed	centroid1	centroid2	radius
# 
# Desired for bigInteract
# chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand


tail -n +3 ${FILE} | awk -v OFS='\t' '{ print "chr"$1, $2, $6, ".", $12, "0", ".", "0", "chr"$1, $2, $3, ".", ".", "chr"$4, $5, $6, ".", "."}' > ${NAME}.bigInteract


tail -n +3 ${FILE} | awk -v OFS='\t' ' $7 != 0 { print "chr"$1, $2, $3, "chr"$4, $5, $6, $12}' > ${NAME}.pairs

# Confirm the rearrangement
echo "Complete"


