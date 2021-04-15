# A script to only extract data from genes mutated in paired samples
# Written by Alex Blain

paired <- BL_mutations %>%
  filter(type=="paired")

paired_frequency <- as.data.frame(table(paired$SYMBOL, paired$Sample_Name))
paired_frequency <- paired_frequency[!(paired_frequency$Freq == 0),]
paired_frequency <- as.data.frame(table(paired_frequency$Var1))
paired_frequency <- paired_frequency[!(paired_frequency$Freq < 3),]

paired_gene_analysis <- BL_mutations[(BL_mutations$SYMBOL %in% paired_frequency$Var1),]

# Mutations per sample
p_count_table1 <- as.data.frame(table(paired_gene_analysis$Sample_Name))

# Mutations per gene per sample
p_count_table2 <-  as.data.frame(table(paired_gene_analysis$SYMBOL, paired_gene_analysis$Sample_Name))
p_count_table2 <- p_count_table2[!(p_count_table2$Freq == 0),]

# Mutations cohort frequency and percentage
p_count_table3 <- as.data.frame(table(p_count_table2$Var1))
p_count_table3$percentage <- (p_count_table3$Freq / 97) * 100
setwd('~/GATK4_analysis2/BL/R_output/')
write_tsv(p_count_table3, "Mutation_BLcohort_pairedanalysis_freq.tsv")

p_count_table3 <- p_count_table3[order(p_count_table3$Freq, decreasing = TRUE), ]
p_count_table3['order'] <- 1:40
order <- p_count_table3 %>%
  select("Var1", "order")



sample_type <- BL_mutations %>%
  select("Sample_Name", "type")
sample_type <- sample_type[!duplicated(sample_type$Sample_Name),]
  
p_count_table4 <-  inner_join(p_count_table2, sample_type, by=c("Var2"="Sample_Name"))
p_count_table5 <- as.data.frame(table(p_count_table4$Var1, p_count_table4$type))
colnames(p_count_table5) <- c("Gene", "Type", "Freq")
p_count_table5 <-  inner_join(p_count_table5, order, by=c("Gene"="Var1"))


plot_p1 <- p_count_table5 %>%
    filter(order <=50, 'Freq' > 2) %>%
  ggplot(p_count_table5, mapping =aes(fill=Type, y=Freq, x=reorder(Gene, order))) +
  geom_bar(position = position_stack(reverse = TRUE), stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Gene") +
  ylab("Number of Samples with Mutation") +
  ggtitle("Number of mutations per sample type")

plot_p1
