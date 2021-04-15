# Script for downstream mutation filtering and visualisation
# Base script written by Dr Masood Zaka and adapted by Alex Blain≈

library(tidyverse)
library(ggplot2)
library(maftools)
library(RColorBrewer)
library(gridExtra)
library(ComplexHeatmap)
# MT2
# vcf

setwd(path_to_paired_vcf)

# reads .tsv files format: chr    pos   ref-allele alt_allele  ref_reads       var_reads   total_depth allele_freq


path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2 <- NULL
for(i in 1:length(file.names)){
  data <- read.table(file.names[i], sep = "\t", header = TRUE,as.is = TRUE)
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[3]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(data))
  data <- cbind(Sample_Name, data)
  data$Location <- paste0(data$chr,":",data$pos)
  data <- data[,c("Sample_Name", "Location", "chr","pos", "ref_allele", "alt_allele", "ref_reads", "var_reads", "total_depth", "allel_freq")]
  MT2 <- rbind(MT2, data)
}

dim(MT2)
head(MT2)
table(MT2$Sample_Name)

# MT2
# vep

setwd(path_to_paired_vep)

path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2_veps <- NULL
for (i in 1:length(file.names)){
  data <- read_tsv(file.names[i], col_names = TRUE, col_types = cols(gnomAD_NFE_AF= col_double(),Location=col_character()))
  filtered.data <- data %>%
    select("Location", "Allele", "Consequence",	"Protein_position",	"Amino_acids","Codons","STRAND",	"VARIANT_CLASS",	"SYMBOL","BIOTYPE","SIFT","PolyPhen", "gnomAD_NFE_AF") %>%
    mutate(gnomAD_NFE_AF= replace_na(gnomAD_NFE_AF,0), Consequence = gsub("[,].*","",Consequence)) %>%
    filter(BIOTYPE=="protein_coding", Amino_acids !="-", str_detect(Consequence, "synonymous",negate = TRUE), gnomAD_NFE_AF < 0.01)
  Amino_acid <- str_split_fixed(filtered.data$Amino_acids,pattern = "/",2)
  Protein_position <- str_split_fixed(filtered.data$Protein_position,"/",2)
  AA <- paste0(Amino_acid[,1],Protein_position[,1],Amino_acid[,2])
  filtered.data <- mutate(filtered.data,Amino_acids=AA)
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[3]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(filtered.data))
  filtered.data <- cbind(Sample_Name, filtered.data)
  MT2_veps <- rbind(MT2_veps, filtered.data)
}

head(MT2_veps)
dim(MT2_veps)

####################################################################################################################
# tumour_only
# vcf
####################################################################################################################
setwd(path_to_TO_vcf)

# reads .tsv files format: chr    pos   ref-allele alt_allele  ref_reads       var_reads   total_depth allele_freq

path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2_tumour_only <- NULL
for(i in 1:length(file.names)){
  data <- read.table(file.names[i], sep = "\t", header = TRUE, as.is = TRUE)
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[1]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(data))
  data <- cbind(Sample_Name, data)
  data$Location <- paste0(data$chr,":",data$pos)
  data <- data[,c("Sample_Name", "Location", "chr","pos", "ref_allele", "alt_allele", "ref_reads", "var_reads", "total_depth", "allel_freq")]
  MT2_tumour_only <- rbind(MT2_tumour_only, data)
}

dim(MT2_tumour_only)
head(MT2_tumour_only)

setwd(path_to_TO_vep)

path = "."
file.names <- dir(path, pattern =".tsv")
file.names
MT2_tumour_only_veps <- NULL
for (i in 1:length(file.names)){
  data <- read_tsv(file.names[i], col_types = cols(gnomAD_NFE_AF= col_double(), Location=col_character()))
  colnames(data)
  filtered.data <- data %>%
    select("Location", "Allele", "Consequence",	"Protein_position",	"Amino_acids","Codons","STRAND",	"VARIANT_CLASS",	"SYMBOL","BIOTYPE","SIFT","PolyPhen", "gnomAD_NFE_AF") %>%
    mutate(gnomAD_NFE_AF= replace_na(gnomAD_NFE_AF,0), Consequence = gsub("[,].*","",Consequence)) %>%
    filter(BIOTYPE=="protein_coding", Amino_acids !="-", str_detect(Consequence, "synonymous",negate = TRUE), str_detect(VARIANT_CLASS, "substitution", negate = TRUE), gnomAD_NFE_AF < 0.01)
  Amino_acid <- str_split_fixed(filtered.data$Amino_acids,pattern = "/",2)
  Protein_position <- str_split_fixed(filtered.data$Protein_position,"/",2)
  AA <- paste0(Amino_acid[,1],Protein_position[,1],Amino_acid[,2])
  filtered.data <- mutate(filtered.data,Amino_acids=AA)
  sampleName <- unlist(strsplit(file.names[i], ".", fixed = TRUE))[1]
  print(paste0("* Running analysis for the :" , sampleName))
  Sample_Name  <- rep(sampleName,nrow(filtered.data))
  filtered.data <- cbind(Sample_Name, filtered.data)
  MT2_tumour_only_veps <- rbind(MT2_tumour_only_veps, filtered.data)
}

head(MT2_tumour_only_veps)
dim(MT2_tumour_only_veps)
table(MT2_tumour_only_veps$Sample_Name)

####################################################################################################################
# vcf and vep joining
####################################################################################################################

# MT2 files

MT2 #vcf
MT2_veps #veps

dim(MT2)
dim(MT2_veps)
head(MT2)
head(MT2_veps)

MT2_veps$Location <- gsub("[-].*","",MT2_veps$Location) # trim after - to include the inframe or framshit variants

MuTect2 <- inner_join(MT2_veps, MT2, by = c("Sample_Name"="Sample_Name","Location"="Location"))
dim(MuTect2)
head(MuTect2)

# MT2_tumour_only files

MT2_tumour_only #vcf
MT2_tumour_only_veps #veps

dim(MT2_tumour_only)
dim(MT2_tumour_only_veps)
head(MT2_tumour_only)
head(MT2_tumour_only_veps)

MT2_tumour_only_veps$Location <- gsub("[-].*","",MT2_tumour_only_veps$Location) # trim after - to include the inframe or framshit variants
MuTect2_Tumour_Only <- inner_join(MT2_tumour_only_veps, MT2_tumour_only, by = c("Sample_Name"="Sample_Name","Location"="Location"))
dim(MuTect2_Tumour_Only)
head(MuTect2_Tumour_Only)

# check
table(colnames(MuTect2)==colnames(MuTect2_Tumour_Only)) # true for correct order

mutations <- rbind(MuTect2,MuTect2_Tumour_Only)
table(order(mutations$SYMBOL))

# Further processing

mutations <- mutations %>%
  mutate(Codons=str_replace_all(Codons,"[:lower:]","")) # remove lower letter from the Codon column
mutations <- mutations %>%
  separate(Codons,c("Ref_allele","Alt_allele"),"/")

mutations <- mutations[mutations$var_reads >= 4 & mutations$total_depth >= 10,] #

dim(mutations)
head(mutations)


setwd(path_to_R_output)
write_tsv(mutations, "mutations.tsv")




########################################################################
#### plots
########################################################################

p1<- mutations %>%
  select(Sample_Name) %>%
  group_by(Sample_Name) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(Sample_Name, -n),y=n)) +
  geom_bar(stat = "identity")+
  xlab("Sample Names") +
  ylab("Number of Mutations") +
  ggtitle("Mutations per Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(size=15))

p2<- mutations %>%
  select(chr) %>%
  group_by(chr) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(chr, -n),y=n)) +
  geom_bar(stat = "identity")+
  xlab("Sample Names") +
  ylab("Number of Mutations") +
  ggtitle("Mutations Load by Chromosome") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p3<- mutations %>%
  select(SYMBOL) %>%
  group_by(SYMBOL) %>%
  count() %>%
  filter(n>20) %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(SYMBOL, -n),y=n)) +
  geom_bar(stat = "identity")+
  xlab("Genes Names") +
  ylab("Number of Mutations > 20") +
  ggtitle("Mutations Load by Gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p4<- mutations %>%
  select(Consequence) %>%
  group_by(Consequence) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(Consequence, n),y=n, fill=Consequence)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  coord_flip() +
  ylab("Variant Counts") +
  ggtitle("Variant Classification") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),axis.title.y = element_blank())

p5<- mutations %>%
  select(VARIANT_CLASS) %>%
  group_by(VARIANT_CLASS) %>%
  count() %>%
  arrange(n) %>%
  ggplot(aes(x=reorder(VARIANT_CLASS, n),y=n, fill=VARIANT_CLASS)) +
  geom_bar(stat = "identity",show.legend = FALSE)+
  coord_flip() +
  ylab("Variant Counts") +
  ggtitle("Variant Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),axis.title.y = element_blank())


mycolors <- colorRampPalette(brewer.pal(8,"Dark2"))(31)

p6<- mutations %>%
  select(Sample_Name, allel_freq) %>%
  group_by(Sample_Name) %>%
  ggplot(aes(x=allel_freq, fill=Sample_Name)) +
  geom_density(alpha= 0.5,position = "stack",show.legend = FALSE)+
  scale_fill_manual(values = mycolors) +
  #facet_wrap(~Sample_Name) +
  xlab("Allelic Frequency") +
  ylab("Density") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),)

png("MutationsProfiles.png",units = "in",width = 18,height = 12,res = 300)
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
graphics.off()
