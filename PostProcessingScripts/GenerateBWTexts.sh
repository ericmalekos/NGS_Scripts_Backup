#!/bin/bash

# check number of arguments
if [[ "$#" -lt 1 ]]; then
    echo "Usage: $0 <genome> [prefix]"
    exit 1
fi

# get the genome and prefix from arguments
genome="$1"
prefix="$2"

# Get the absolute path of the current directory
dirpath="$(pwd)"

# if dirpath starts with ~/public_html/ (or the resolved path of that directory), remove that part for the URL
home_public_html="$(realpath ~/public_html)"
urlpath="${dirpath#$home_public_html/}"

# loop through each .bw file in the directory
for bwfile in "$dirpath"/*.bw; do
    # get just the filename without path
    bwfilename="$(basename "$bwfile")"
    
    # construct the output text
    text="track type=bigWig name='${prefix}_${bwfilename}'  description='${prefix}_${bwfilename}'  genome=${genome}  color=0,0,0 altColor=0,0,0  priority=0.034  autoScale=on   viewLimits=1:200  yLineOnOff=off yLineMark=100   visibility=full  maxHeightPixels=100:50:20  bigDataUrl=http://public.gi.ucsc.edu/~emalekos/${urlpath}/${bwfilename}"

    # output the text to a file with the same name but .txt extension in the same directory
    echo "$text" > "${bwfile%.bw}.txt"
done

echo "Done generating text files!"
