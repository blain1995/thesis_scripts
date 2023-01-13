#!/bin/bash 

# Work out the coverage of a specific exon and compare with the GC content of the gene
# If the gene is on the reverse strand, please set the reverse object to 'true', otherwise set to 'false'

# Load required modules
# module load Anaconda3 Python SAMtools
# source activate pythonAB

# reverse='true'

# Run this if the GC content of the gene has not already been calculated
# python ./gene_gc.py

# Run samtools depth to get coverage
samtools depth -a -f coverage_input.txt -H -r 13:41133660-41134997 -o foxo1_exon2_depth.txt

number=$(wc -l coverage_input.txt | awk '{ print $1 }')

echo $number

