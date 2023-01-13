# Alex Blain 2022
# Use to calculate average depth at each specified position

import pandas as pd

# Read in files
df = pd.read_csv("foxo1_exon2_depth.txt", sep="\t")
gene = pd.read_csv("foxo2_gc.tsv", sep="\t")

# Get average and median depth
df['avg_depth'] = df.iloc[:, 2:].mean(axis=1)
df['median_depth'] = df.iloc[:,2:75].median(axis=1)

# Left join with the gc file
gene = gene.merge(df, on='POS')
print(gene.head())

gene.to_csv("avg_depth_exon2.tsv", sep="\t", index=False)
