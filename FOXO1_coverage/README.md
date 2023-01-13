These scripts were used to assess the GC content and sequencing coverage in section 3.2.9 of my PhD thesis.

To run these scripts you first need to download the fasta sequence of the exon of interest from 
[ensembl](https://useast.ensembl.org/) and store as a text file. The text file is then used as input to the 
[gene_gc.py](https://github.com/blain1995/thesis_scripts/FOXO1_coverage/gene_gc.py) script, which builds a python 
dictionary including the 
genomic co-ordinate, base and GC boolean at each position of the exon. 

The next step is to run [coverage.sh](https://github.com/blain1995/thesis_scripts/FOXO1_coverage/coverage.sh) which 
uses samtools depth to determine the sequencing depth for each sample at every given base of the exon. Finally 
[avg_depth.py](https://github.com/blain1995/thesis_scripts/FOXO1_coverage/avg_depth.py) is used to calculate the 
average depth across all samples at every given base.   
