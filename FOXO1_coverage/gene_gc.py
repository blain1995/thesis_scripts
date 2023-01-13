# Alex Blain 2022
# Use to determine the GC content of a sequence
# If gene on reverse strand set reverse to True

import pandas as pd

# UPDATE BEFORE RUNNING
reverse = True

# Read in the sequence text file as a string
file = open("foxo1_exon2.txt")
sequence = file.read().replace("\n", "")
file.close

# Checks
print(sequence)
print(type(sequence))

length = len(sequence)
print(length)

# Make an empty list to append to
gene_gc = []

# Set coordindare to the genomic co-ordinate of the start of the gene
coordinate = 41134997

for base in range(0, length - 1):
	if sequence[base] == 'G' or sequence[base] == 'C':
		GC = 1
	else:
		GC = 0

	sample_dict = {'position': base+1,
	'base':sequence[base],
	'GC_content': GC,
	'POS': coordinate
	}

	gene_gc.append(sample_dict)

	if reverse == True:
		coordinate -= 1
	else:
		coordinate += 1

df = pd.DataFrame(gene_gc)

print(df.head())
df.to_csv("foxo2_gc.tsv", sep="\t", index=False)	

