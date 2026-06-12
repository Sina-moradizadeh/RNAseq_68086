setwd("E:/RNAseq_GSE68086")
dir()

library(tidyverse)

# Data preparation, Meta Data

data <- read.delim("GSE68086_TEP_data_matrix.txt", sep = "\t")
View(data)
colnames(data)

rownames(data) <- data[, 1]
rownames(data)
data <- data[, -1]
ncol(data)


data_sample_filter <- data |>
select(-starts_with("Type.Unknown"))
ncol(data_sample_filter)
write.table(data_sample_filter, "DataSampleFiltere.txt", sep = "\t", quote = F)


samples <- colnames(data_sample_filter)
length(samples)
group <- rep(NA, length(samples))

group[grepl("^HD|^Control", samples)] <- "Health"
group[grepl("Breast|BrCa", samples)] <- "Breast"
group[grepl("CRC", samples)] <- "Colorectal"
group[grepl("GBM", samples)] <- "Brain"
group[grepl("VU383Platelet.hiseq|VU394Platelet.hiseq", samples)] <- "Brain"
group[grepl("Lung|Vumc|MGH", samples)] <- "Lung"
group[grepl("Pancr|Panc", samples)] <- "Pancreas"
group[grepl("Liver", samples)] <- "Liver"
group[grepl("Chol", samples)] <- "Cholangiocarcinoma"


metadata <- data.frame(
  Samples = samples,
  Group = factor(group)
)

table(metadata$Group)
nrow(metadata)
view(metadata)

write.table(metadata, "MetaData.txt", sep = "\t", quote = F)

# Gene Filter(gene with low count)

gene_filter <- data_sample_filter[which(rowSums(data_sample_filter) >= 10), ]

nrow(data_sample_filter)
nrow(gene_filter)

write.table(gene_filter, "gene_filter.txt", sep = "\t", quote = F)

# Annotation / 1

annot <- read_tsv("Human.GRCh38.p13.annot.tsv")

annot_1 <- annot |>
  select(c(2,5,6)) |>
  distinct(EnsemblGeneID, .keep_all = TRUE)
view(annot_1)

write.table(annot_1, "Annot_1", sep = "\t", quote = F)

gene_filter_z <- gene_filter |>
  rownames_to_column("EnsemblGeneID") |>
  left_join(annot_1)

write.table(gene_filter_z, "gene_filter_z.txt", sep = "\t", quote = F)

nrow(gene_filter_z)
nrow(gene_filter)

view(gene_filter)
view(gene_filter_z)
table(gene_filter_z$Symbol)
table(gene_filter_z$GeneType)
sum(is.na(gene_filter_z$Symbol))

gene_filter_z <- read.delim("gene_filter_z.txt", sep = "\t", header = T)
view(gene_filter_z)

unkown_genes <- gene_filter_z |>
  filter(is.na(Symbol))
view(unkown_genes)

# Annotation / 2

gene_filter <- read.delim("gene_filter.txt", sep = "\t", header = T)
dir()
view(gene_filter)
gene_En <- rownames(gene_filter)
nrow(gene_filter)
view(gene_En)

library(org.Hs.eg.db)
library(AnnotationDbi)

symbol <- mapIds(
  org.Hs.eg.db,
  keys = gene_En,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(symbol)
view(symbol)

table(symbol)
length(symbol)
sum(is.na(symbol))

symbol_2 <- data.frame(
  gene_En,
  symbol
)
view(symbol_2)
sum(is.na(symbol_2$symbol))

annot_package <- gene_filter |>
  rownames_to_column(var = "gene_En") |>
  left_join(symbol_2) |>
  relocate(symbol, .after = gene_En)
  
nrow(annot_package)
view(annot_package)
colnames(annot_package)

write.table(annot_package, "Annotation with package.txt", quote = F, sep = "\t")

# Plot
pdf("Box_Plot_befor.pdf", width = 8, height = 6)
boxplot(log2(gene_filter + 1),
        las = 2,
        outline = F
        )
dev.off()
  
pdf("pca_befor.pdf", width = 8, height = 6)
pca <- prcomp(t(log2(gene_filter + 1)))
plot(pca$x[, 1], pca$x[, 2])
dev.off()

# DESeq2

library(DESeq2)
dir()
meta_data <- read.delim("MetaData.txt", sep = "\t")
rownames(meta_data) <- meta_data[, 1]

table(meta_data$Group)
colnames(meta_data)
rownames(meta_data)
view(meta_data)
identical(meta_data$Samples, colnames(gene_filter))

colnames(gene_filter)
ncol(gene_filter)
nrow(meta_data)

d1 <- DESeqDataSetFromMatrix(countData = gene_filter,
                             colData = meta_data,
                             design = ~ Group)
data_vest <- vst(d1, blind = TRUE)
data_vest_dataframe <- data.frame(assay(data_vest))
view(data_vest_dataframe)
write.table(data_vest_dataframe, "vest data .txt", quote = F, sep = "\t")

pdf("Box_Plot_After.pdf", width = 8, height = 6)
boxplot(data_vest_dataframe,
        las = 2,
        outline = F
)
dev.off()

pdf("pca_Plot_after.pdf", width = 8, height = 6)
pca <- prcomp(t(data_vest_dataframe))
plot(pca$x[, 1], pca$x[, 2])
dev.off()

d2 <- estimateSizeFactors(d1)

d3 <- DESeq(d2)

resultsNames(d3)


res_breast <- results(
  d3,
  contrast = c(
    "Group",
    "Breast",
    "Health"
  )
)
write.table(res_breast, "BreastVsHealth.txt", sep = "\t")


res_brain <- results(
  d3,
  contrast = c(
    "Group",
    "Brain",
    "Health"
  )
)
write.table(res_brain, "BrainVsHealth.txt", sep = "\t")


res_col <- results(
  d3,
  contrast = c(
    "Group",
    "Colorectal",
    "Health"
  )
)
write.table(res_col, "ColorectalVsHealth.txt", sep = "\t")


res_lung <- results(
  d3,
  contrast = c(
    "Group",
    "Lung",
    "Health"
  )
)
write.table(res_lung, "LungVsHealth.txt", sep = "\t")


res_pancreas <- results(
  d3,
  contrast = c(
    "Group",
    "Pancreas",
    "Health"
  )
)
write.table(res_pancreas, "PancreasVsHealth.txt", sep = "\t")

pdf("Plot pca end.pdf", width = 8, height = 6)
plotPCA(
  data_vest,
  intgroup = "Group"
)
dev.off()














