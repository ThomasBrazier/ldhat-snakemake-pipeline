#!/usr/bin/env python
"""
Script to cut large list in concatenated list elements with
        overlapping options

Usage:
    ./<script.py> inputfile outputfile size overlap

Where:
    inputfile: input file name
    outputfile: path to output file with prefix
    size:     is the number of SNPs to keep
    overlap:  is the number of overlapping SNPs

For example, 2,000 total SNPs, with 200 SNPs at the beginning/end shared with previous/next chunk
"""

import sys
import gzip

try:
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
    size = int(sys.argv[3])
    overlap = int(sys.argv[4])
except:
    print(__doc__)
    sys.exit(1)

def get_overlapped_chunks(text, chunksize, overlapsize):  
    return [ lines[a:a+chunksize] for a in range(0,len(lines), chunksize-overlapsize)]

lines = open(inputfile, 'r').readlines() 

chunks = get_overlapped_chunks(lines, size, overlap)

length = len(chunks) 
filenames=[f"{outputfile}/batch_{x+1}.pos" for x in range(length) ]

open(f"{outputfile}/nbatch_split","w").write(str(length))

for chunk, files in zip(chunks, filenames):
    with open(files, 'w') as output:
        output.write("".join(chunk))


