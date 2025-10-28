# Data import 
counts <- read.table("all_samples_counts.csv",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)


suppressPackageStartupMessages(library(DESeq2)) # to load DESeq2 and suppress the long startup message

# Make metadata
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "treated", "treated")
)

# Convert condition to factor
coldata$condition <- factor(coldata$condition)

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
print(dds)


