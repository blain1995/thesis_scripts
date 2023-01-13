# Differential expression analysis detailed in chapter six of my thesis 
# Alex Blain 2021

library(devtools)
install_github("mikelove/tximport")
library(tximport)
library(tidyverse)
library(DESeq2)
library(GenomicFeatures)
library(ggsignif)
library(ComplexHeatmap)
library(magrittr)

#### Load sample info and generate metadata file ####
setwd("~/Documents/RNA-seq/")

# Format input data
sample_table <- read_tsv("sample-table-md2.txt")
sample_table <- mutate(sample_table, 
                       files=paste0("~/Documents/RNA-seq/counts_gcbias_unsorted/", ID, "/quant.sf")) %>%
  mutate(., condition=paste0(diagnosis, "_", relapse, "_", sample_type))

files = pull(sample_table, files)
names(files) <- sample_table$names
file.exists(files)

coldata <- data.frame(files, names=sample_table$names, 
                      condition=sample_table$condition, stringsAsFactors = FALSE)


#### Make the transcript 2 gene file ####
txdb <- makeTxDbFromGFF(file="gencode.v37lift37.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")

columns(txdb)
# Should return a 1:1 mapping
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
# write.table(tx2gene, file="tx2gene.txt", sep="\t")

test <- tx2gene
test2 <- test
test2$Transcript.stable.ID <- str_split(test2$TXNAME, "\\.", simplify=T)[,1]
test2$Gene.stable.ID <- str_split(test2$GENEID, "\\.", simplify=T)[,1]

mart <- read.table("mart_export.txt", sep="\t", header=TRUE)
genes <- mart %>% dplyr::select(Gene.stable.ID.version, Gene.name)
# write.table(genes, file="genes.txt", sep = "\t")

annotation <- left_join(test2, mart)
# write.table(annotation, file="annotation.txt", sep="\t")

genes2 <- annotation %>% dplyr::select(., GENEID, Gene.name) %>% unique()
colnames(genes2) <- c("gene", "gene_name")
# write.table(genes2, file="genes2.txt", sep="\t")

#### Running DESeq2 ####
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData=coldata,
                                   design=~condition)
dds <- DESeq(ddsTxi)
# write.table(counts(dds), file="original_counts.txt")
resultsNames(dds)

counts <- counts(dds)
normalised <- counts(dds, normalized=TRUE)
# write.table(normalised, file="DESeq2_normalised_counts.txt", sep="\t")

rlog_out <- rlog(dds, blind=FALSE)
rlog_counts <- assay(rlog_out)
vsd <- varianceStabilizingTransformation(dds)
# write.table(assay(vsd), file="varianceStabilisingTransformedCounts.txt", sep="\t")

#### BL Rel vs no Results ####
# dir.create("R_output_final")
setwd("~/Documents/RNA-seq/R_output_final/")
res <- results(dds, contrast = c("condition", "BL_No_Diagnosis", "BL_Yes_Diagnosis"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

# write.table(res_df, file="BLNoDiagnosis.vs.BLYesDiagnosis.txt", sep="\t")
BL_res_ann <- dplyr::left_join(res_df, genes2, by=c("gene")) %>% na.omit()
# write.table(BL_res_ann, file="BLNoDiagnosis.vs.BLYesDiagnosis_Annotated.txt", sep="\t")

res0.5 <- results(dds,contrast = c("condition", "BL_No_Diagnosis", "BL_Yes_Diagnosis"), alpha=0.05)
sum(res$padj < 0.05, na.rm=TRUE)
res0.5_df <- as.data.frame(res0.5) %>% filter(padj < 0.05)
res0.5_df$gene <- rownames(res0.5_df)
# write.table(res0.5, file="BLNoDiagnosis.vs.BLYesDiagnosis_padj0.5.txt", sep="\t")

BL_sig_ann <- dplyr::left_join(res0.5_df, genes2, by=c("gene")) %>% na.omit()
# write.table(BL_res_ann, file="BLNoDiagnosis.vs.BLYesDiagnosis_Annotated.txt", sep="\t")

BL_top20 <- BL_sig_ann[order(BL_sig_ann$padj),][1:20,]

#### BL Rel diag vs rel Results ####
resREL <- results(dds, contrast = c("condition", "BL_Yes_Diagnosis", "BL_Yes_Relapse"))
resREL_df <- as.data.frame(resREL)
resREL_df$gene <- rownames(resREL_df)
summary(resREL)
sum(resREL$padj < 0.1, na.rm=TRUE)

# write.table(resREL_df, file="BLYesDiagnosis.vs.BLYesRelapse.txt", sep="\t")
BL_resREL_ann <- dplyr::left_join(resREL_df, genes2, by=c("gene")) %>% na.omit()
# write.table(BL_resREL_ann, file="BLYesDiagnosis.vs.BLYesRelapse.txt", sep="\t")

resREL0.5 <- results(dds,contrast = c("condition", "BL_Yes_Diagnosis", "BL_Yes_Relapse"), alpha=0.05)
sum(resREL$padj < 0.05, na.rm=TRUE)
resREL0.5_df <- as.data.frame(resREL0.5) %>% filter(padj < 0.05)
resREL0.5_df$gene <- rownames(resREL0.5_df)
# write.table(resREL0.5, file="BLYesDiagnosis.vs.BLYesRelapse_padj0.5.txt", sep="\t")

BL_RELsig_ann <- dplyr::left_join(resREL0.5_df, genes2, by=c("gene")) %>% na.omit()
# write.table(BL_res_ann, file="BLYesDiagnosis.vs.BLYesRelapse_Annotated.txt", sep="\t")

#### BLL-11q vs BL Results ####
res11q <- results(dds, contrast = c("condition", "BLL.11q_No_Diagnosis", "BL_No_Diagnosis"))
summary(res11q)
sum(res11q$padj < 0.1, na.rm=TRUE)
res11q_df <- as.data.frame(res11q)
res11q_df$gene <- rownames(res11q_df)

# write.table(res11q_df, file="BLL11q.vs.BL.txt", sep="\t")
BLL11q_res_ann <- dplyr::left_join(res11q_df, genes2, by=c("gene")) %>% na.omit()
# write.table(BLL11q_res_ann, file="BLL11q.vs.BL_Annotated.txt", sep="\t")

sig11q <- results(dds, contrast = c("condition", "BLL.11q_No_Diagnosis", "BL_No_Diagnosis"), alpha=0.05)
BLL11q0.5 <- as.data.frame(sig11q) %>% filter(padj < 0.05)
BLL11q0.5$gene <- rownames(BLL11q0.5)
# write.table(BLL11q0.5, file="BLL11q.vs.BL_padj0.5.txt", sep="\t")

BLL11q_sig_ann <- dplyr::left_join(BLL11q0.5, genes2, by=c("gene")) %>% na.omit()
# write.table(BLL11q_sig_ann, file="BLL11q.vs.BL_padj0.5_Annotated.txt", sep="\t")

BLL11q_annotated <- dplyr::left_join(BLL11q0.5, genes2, by=c("gene")) %>% na.omit()

# BLL-11q heatmap
heatmap_genes <- c("CREB5", "ETS1", "IL10RA", "CD83", "PRKCD", "IFI35", "TFEB",
                   "STAT3", "BCL2L1", "IFIT1", "IRF7", "BCL6", "OAS3", "FEZ1",
                  "LAMP3", "BACH2", "NFRKB", "MYD88", "ID3", "MYC")

heatmap_genes <- c("IL10RA", "STAT3", "STAT5A", "TNF", "BLVRA", "IL10", "JAK1",
                   "STAT1")

c("HMOX1", "BLVRB", "IL10RB", "IL1A", "IL6")

ensembl_id <- annotation %>% filter(Gene.name %in% heatmap_genes) %>%
  dplyr::select(GENEID, Gene.name) %>% unique()
rownames(ensembl_id) <- ensembl_id$GENEID

heatmap_counts <- rlog_counts[rownames(rlog_counts) %in% ensembl_id$GENEID, 
                             colnames(rlog_counts) %in% BLL11q_GSEA$names] 

annotated_heatmap <- merge(ensembl_id, heatmap_counts, by=0)[,-(1:2)] %>%
  set_rownames(.$Gene.name) %>% dplyr::select(-Gene.name) %>%
  as.matrix(.)


heatmap_order <- read.table("heatmap_order.txt", sep="\t")

heatmap.scaled <- t(apply(annotated_heatmap, 1, scale)) %>%
  set_colnames(., colnames(annotated_heatmap))
  
heatmap.scaled.ordered <- heatmap.scaled[match(heatmap_genes, rownames(heatmap.scaled)), match(heatmap_order$V1, colnames(heatmap.scaled))]

ha <- HeatmapAnnotation(Diagnosis=rep(c("BLL-11q", "BL"), times=c(3, 40)),
                                      col=list(Diagnosis=c("BLL-11q"="#67AD57",
                                                           "BL"="#4A7CB3")))
Heatmap(heatmap.scaled.ordered, name="z-score", cluster_columns=FALSE, cluster_rows=FALSE,top_annotation = ha)

vsd_counts <- assay(vsd)

df_test <- as.data.frame(colData(dds))
#### BLL-11q vs DLBCL Results ####
res11qDLBCL <- results(dds, contrast = c("condition", "BLL.11q_No_Diagnosis", "DLBCL_No_Diagnosis"))
summary(res11qDLBCL)
sum(res11qDLBCL$padj < 0.1, na.rm=TRUE)
res11qDLBCL_df <- as.data.frame(res11qDLBCL)
res11qDLBCL_df$gene <- rownames(res11qDLBCL_df)

# write.table(res11qDLBCL_df, file="BLL11q.vs.DLBCL.txt", sep="\t")
BLL11qDLBCL_res_ann <- dplyr::left_join(res11qDLBCL_df, genes2, by=c("gene")) %>% na.omit()
# write.table(BLL11qDLBCL_res_ann, file="BLL11q.vs.DLBCL_Annotated.txt", sep="\t")

sig11qDLBCL <- results(dds, contrast = c("condition", "BLL.11q_No_Diagnosis", "DLBCL_No_Diagnosis"), alpha=0.05)
BLL11qDLBCL0.5 <- as.data.frame(sig11qDLBCL) %>% filter(padj < 0.05)
BLL11qDLBCL0.5$gene <- rownames(BLL11qDLBCL0.5)
# write.table(BLL11qDLBCL0.5, file="BLL11q.vs.DLBCL_padj0.5.txt", sep="\t")

BLL11qDLBCL_sig_ann <- dplyr::left_join(BLL11qDLBCL0.5, genes2, by=c("gene")) %>% na.omit()
# write.table(BLL11qDLBCL_sig_ann, file="BLL11q.vs.DLBCL_padj0.5_Annotated.txt", sep="\t")

BLL11qDLBCL_annotated <- dplyr::left_join(BLL11qDLBCL0.5, genes2, by=c("gene")) %>% na.omit()
BLL11q_BL_DLBCL <- coldata %>% filter(condition == "BL_No_Diagnosis" | condition == "BLL-11q_No_Diagnosis" | condition == "DLBCL_No_Diagnosis")
list <- BLL11q_BL_DLBCL$names
subtype_counts <- counts[, list]

sg <- c("CKLF", "CCND3", "CREB5", "CD53", "IL16", "ETS1", "TIGIT")
sg_df <- BLL11qDLBCL_res_ann[BLL11q_res_ann$gene_name %in% sg, ]

plotCounts(dds, gene="ENSG00000077238.14_10")

subtype_counts <- normalised[, list]
subtype_counts$gene <- rownames(subtype_counts)
subtype_annotated <- dplyr::left_join(BLL11q0.5, genes2, by=c("gene")) %>% na.omit()

#### GSEA set up ####
coldata
BL_GSEA <- coldata %>% filter(condition == "BL_No_Diagnosis" | condition == "BL_Yes_Diagnosis")
BL_GSEA_counts <- normalised[,BL_GSEA$names]
BL_GSEA_counts <- as.data.frame(BL_GSEA_counts)
BL_GSEA_counts$gene <- rownames(BL_GSEA_counts)
BL_GSEA_counts_ann <- dplyr::left_join(BL_GSEA_counts, genes2, by=c("gene")) %>% na.omit()
# write.table(BL_GSEA_counts_ann, file="BL_GSEA.txt", sep="\t")
BL_GSEA_description <- BL_GSEA %>% dplyr::select(names, condition)
# write.table(BL_GSEA_description, file="BL_GSEA_description.txt", sep="\t", row.names = FALSE)

BLL11q_GSEA <- coldata %>% filter(condition == "BLL-11q_No_Diagnosis" | condition == "BL_No_Diagnosis")
BLL11q_GSEA_counts <- normalised[,BLL11q_GSEA$names]
BLL11q_GSEA_counts <- as.data.frame(BLL11q_GSEA_counts)
BLL11q_GSEA_counts$gene <- rownames(BLL11q_GSEA_counts)
BLL11q_GSEA_counts_ann <- dplyr::left_join(BLL11q_GSEA_counts, genes2, by=c("gene")) %>% na.omit()
# write.table(BLL11q_GSEA_counts_ann, file="BLL11q_GSEA.txt", sep="\t", row.names=FALSE)
BLL11q_GSEA_description <- BLL11q_GSEA %>% dplyr::select(names, condition)
# write.table(BLL11q_GSEA_description, file="BLL11q_GSEA_description.txt", sep="\t", row.names = FALSE)

#### Get counts and annotations for a specific gene ####
log2counts <- log2(counts +1)
subset_counts <- as.data.frame(log2counts[,BL_GSEA$names]) %>% 
  rownames_to_column(.,) %>% left_join(., genes2, by=c("rowname" = "gene")) %>%
  na.omit() %>% dplyr::select(-rowname)


TFRC <- as.data.frame(t(subset_counts %>% filter(gene_name == "TFRC"))[1:46,]) %>%
  rownames_to_column(.,) %>% magrittr::set_colnames(., c("names", "counts")) %>%
  left_join(., BL_GSEA, by="names") %>% dplyr::select(., -files) 

TFRC$counts <- as.numeric(TFRC$counts)

ggplot(TFRC, aes(x=condition, y=counts, fill=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  geom_signif(comparisons = list(c("BL_No_Diagnosis", "BL_Yes_Diagnosis")),
              map_signif_level = TRUE,
              tip_length = 0) +
  theme_bw() +
  theme(legend.position = "none") + 
              #axis.text=element_text(size=14),
              #axis.title=element_text(size=18)) +
          ylab("RNA-seq Log2 counts") +
          ylim(10,20) +
          scale_fill_manual(values=c("#377eb8", "#e51b18"))

# Compare to NanoString data if needed
NanoString <- read_csv("NanoString_log2counts.csv")
NanoString_annotations <- read_csv("annotations.csv")

NanoString_coldata = NanoString_annotations %>% dplyr::select(., `Sample ID`, condition)
NanoString_subset <- NanoString_coldata %>% filter(condition == "BL_No_Diagnosis" | condition == "BLL-11q_No_Diagnosis" | condition == "DLBCL_No_Diagnosis")

TFRC_NanoString <- NanoString %>% dplyr::select("Sample", "CCND3-mRNA") %>% 
  filter(Sample %in% NanoString_subset$`Sample ID`) %>%
  left_join(., NanoString_subset, by=c("Sample"="Sample ID")) %>%
  magrittr::set_colnames(., c("names", "counts", "condition"))

ggplot(TFRC_NanoString, aes(x=condition, y=counts, fill=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  geom_signif(comparisons = list(c("BL_No_Diagnosis", "BLL-11q_No_Diagnosis"), c("BLL-11q_No_Diagnosis","DLBCL_No_Diagnosis")),
              map_signif_level = TRUE)+
              # tip_length = 0) +
  theme_bw() +
  theme(legend.position = "none", 
  axis.text=element_text(size=14),
  axis.title=element_text(size=16)) +
  ylab("NanoString CCND3 Log2 counts") +
  ylim(6.5,11.5) +
  scale_x_discrete(labels=c("BL", "BLL-11q", "DLBCL")) +
  scale_fill_manual(values=c("#377eb8", "#e51b18", "#67AD57"))

TFRC$platform <- "RNA_seq"
TFRC_NanoString$platform <- "NanoString"

TFRC_all_counts <- rbind(TFRC, TFRC_NanoString) 
ggplot(TFRC_all_counts, aes(x=condition, y=counts, fill=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  geom_signif(comparisons = list(c("BL_No_Diagnosis", "BL_Yes_Diagnosis")),
              map_signif_level = TRUE,
              tip_length = 0) +
  facet_wrap(~platform) +
  theme_bw() +
  theme(legend.position = "none") + 
  #axis.text=element_text(size=14),
  #axis.title=element_text(size=18)) +
  ylab("log2 counts") +
  scale_fill_manual(values=c("#377eb8", "#e51b18"))


#### 14/486 removal check ####
files2 <- files[files !="~/Documents/RNA-seq/counts_gcbias_unsorted/68_14_486_Vikki_Rand_S68/quant.sf"] 
coldata2 <- coldata %>% filter(names != "14_486")
txi2 <- tximport(files2, type="salmon", tx2gene=tx2gene)
ddsTxi2 <- DESeqDataSetFromTximport(txi2,
                                   colData=coldata2,
                                   design=~condition)

dds2 <- DESeq(ddsTxi2)
res2 <- results(dds2, contrast = c("condition", "BL_No_Diagnosis", "BL_Yes_Diagnosis"))
res_df2 <- as.data.frame(res2)
res_df2$gene <- rownames(res_df2)
res2_ann <- dplyr::left_join(res_df2, genes2, by=c("gene")) %>% na.omit()


#### IL10 BL investigation ####
samples <- c("15_476", "16_817", "18_622", "14_228", "20_3044")
BLw11q <- c("4_911", "15_338", "22_9")

IL10_samples <- sample_table %>%
  mutate(IL10 = case_when(names %in% samples ~ 'IL10',
                          names %in% BLw11q ~ 'w11q',
                          !(names %in% samples | names %in% BLw11q) ~ 'no'),
         IL10_condition = paste0(diagnosis, "_", IL10))

IL10_comp <- IL10_samples %>% filter(., IL10_condition == "BL_IL10" | IL10_condition == "BL_no" |
                                       IL10_condition == "BLL-11q_no" | IL10_condition == "BL_w11q") %>%
  dplyr::select(names, IL10_condition)
  

log2counts <- log2(counts +1)
subset_counts <- as.data.frame(log2counts[,IL10_comp$names]) %>% 
  rownames_to_column(.,) %>% left_join(., genes2, by=c("rowname" = "gene")) %>%
  na.omit() %>% dplyr::select(-rowname)


IL10RA <- as.data.frame(t(subset_counts %>% filter(gene_name == "MYC"))[1:53,]) %>%
  rownames_to_column(.,) %>% magrittr::set_colnames(., c("names", "counts")) %>%
  left_join(., IL10_comp, by="names") 

IL10RA$counts <- as.numeric(IL10RA$counts)

ggplot(IL10RA, aes(x=IL10_condition, y=counts, fill=IL10_condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  # geom_signif(comparisons = list(c("BL_No_Diagnosis", "BL_Yes_Diagnosis")),
              # map_signif_level = TRUE,
              # tip_length = 0) +
  theme_bw() +
  theme(legend.position = "none") + 
  #axis.text=element_text(size=14),
  #axis.title=element_text(size=18)) +
  ylab("RNA-seq Log2 counts") +
  # ylim(10,20) +
  scale_fill_manual(values=c("#377eb8", "#e51b18", "green", "orange"))

