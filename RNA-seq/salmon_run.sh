#!/bin/bash
reads=/home/fastq
index=/home/salmon_index
counts=/home/counts

for a in $(cat sample_list.txt)
	do
	salmon quant -i ${index} --libType A --validateMappings \
	-1 <(gunzip -c $reads/${a}_L001_R1_001.fastq.gz $reads/${a}_L002_R1_001.fastq.gz $reads/${a}_L003_R1_001.fastq.gz $reads/${a}_L004_R1_001.fastq.gz) \
	-2 <(gunzip -c $reads/${a}_L001_R2_001.fastq.gz $reads/${a}_L002_R2_001.fastq.gz $reads/${a}_L003_R2_001.fastq.gz $reads/${a}_L004_R2_001.fastq.gz) \
	-p 16 \
	-o ${counts}/${a}
	done
