library(Rsubread)

# Run featureCounts on all 4 SAM files together
out_FeatCount <- featureCounts(
  files = c("conA.sam", "conA1.sam", "conB1.sam", "conB2.sam"),
  annot.ext = "Copy of GCF_000146045.2_R64_genomic.gtf",
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  nthreads = 8  # adjust based on your CPU
)

# Extract counts matrix
counts_matrix <- as.data.frame(out_FeatCount$counts)

# Add gene_id column from annotation
counts_matrix <- cbind(gene_id = rownames(counts_matrix), counts_matrix)

# Rename sample columns (optional, for clarity)
colnames(counts_matrix) <- c("gene_id", "conA", "conA1", "conB1", "conB2")

# Save combined table to file
write.table(
  counts_matrix,
  file = "all_samples_counts.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# View first rows
head(counts_matrix)

rownames(counts_matrix)<-NULL

LibSize_WT1 = sum(counts_matrix$conA)
LibSize_WT2 = sum(counts_matrix$conA1)
LibSize_MT1=sum(counts_matrix$conB1)
LibSize_MT2=sum(counts_matrix$conB2)

sf_WT1 = LibSize_WT1/1000000
sf_WT2 = LibSize_WT2/1000000
sf_MT1=LibSize_MT1/1000000
sf_MT2=LibSize_MT2/1000000

counts_matrix$WT1_rpm = counts_matrix$conA/ sf_WT1
counts_matrix$WT2_rpm = counts_matrix$conA1 / sf_WT2

counts_matrix$MT1_rpm = counts_matrix$conB1/ sf_MT1
counts_matrix$MT2_rpm = counts_matrix$conB2 / sf_MT2

df = read.table("gene_length.txt",header =TRUE)

gene_lgth_kb = (df$Length.bp./1000)

RPKM_WT1 = counts_matrix$WT1_rpm/gene_lgth_kb
RPKM_WT2 = counts_matrix$WT2_rpm/gene_lgth_kb
RPKM_MT1=counts_matrix$MT1_rpm/gene_lgth_kb
RPKM_MT2=counts_matrix$MT2_rpm/gene_lgth_kb

length(counts_matrix$gene_id)
length(RPKM_WT1)
length(RPKM_WT2)
length(RPKM_MT1)
length(RPKM_MT2)

print(counts_matrix$gene_id)

print(RPKM_MT1[6459:6465])

count_matrix_unique <- count_matrix[!duplicated(count_matrix$gene_id), ]

