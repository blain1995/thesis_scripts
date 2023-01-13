# RNA-Seq analysis using Salmon and DESeq2

## Analysis set up

There are many methods developed for the mapping and quantification of RNA sequencing data. For my analysis I chose to 
use the Salmon package developed by the [COMBINE-lab](https://github.com/COMBINE-lab/salmon). For my analysis I decided 
to run this using the combine lab's [docker image](https://hub.docker.com/r/combinelab/salmon).

To run the analysis you first need to obtain the reference files from the genome build of interest. For this project I 
used [Gencode](https://www.gencodegenes.org). This analysis requires a transcriptome file used for pseudoalignment and 
a genome file used as a decoy. For the analysis carried out as part of chapter six in my thesis I used the [v37 GRCh37 
liftover](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/) as my reference. The 
transcriptome file is named "gencode.v37lift37.transcripts.fa.gz" and the genome build is 
"GRCh37.primary_assembly.genome.fa.gz".

---

## Building the reference index

The first step is to build the dedcoy-aware reference transcriptome. This can be done by using the code below:

```bash
# Create the decoy file
grep "^>" <(gunzip -c GRCh37.primary_assembly.genome.fa.gz) | \
cut -d "	" -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# Concatenate genome/transcriptome reference
cat gencode.v37lift37.fa GRCh37.primary_assembly.genome.fa > gentrome.fa

# Build the decoy-aware salmon index
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index--gencode
```

This reference index is then used in subsequent salmon mapping experiments.

## Analysis

After using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to confirm the quality of each fastq 
file is sufficient to continue, fastq files can be mapped using the 
[salmon_run.sh](https://github.com/blain1995/thesis_scripts/tree/main/RNA-seq/salmon_run.sh). The output is quant.sf 
files for each sample. Downstream analysis of quant.sf files is carried out using 
[R_analysis_final.R](https://github.com/blain1995/thesis_scripts/tree/main/RNA-seq/R_analysis_final.R), which is used 
to make the [transcript to gene](https://github.com/blain1995/thesis_scripts/tree/main/RNA-seq/tx2gene.txt) file and 
to carry out differential expression analysis.


 
