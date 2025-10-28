# Install if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures"))

library(rtracklayer)

# Import the GTF
gtf <- import("genome_fixed.gtf")

# Extract relevant info
gtf_df <- as.data.frame(gtf)

# Look at available columns
head(colnames(gtf_df))

library(dplyr)

gene_length_df <- gtf_df %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  select(gene_id, start, end) %>%
  mutate(length_bp = end - start + 1)

print(gene_length_df)

counts <- read.table("all_samples_counts.csv",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)
# If your counts_df has gene_id as rownames, first turn it into a column:
counts <- counts %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "gene_id")



merged_df <- counts %>%
  inner_join(gene_length_df, by = "gene_id")


# Now merged_df has counts + length_bp
head(merged_df)

# Extract counts matrix
counts <- as.matrix(merged_df[,2:5])
gene_length_kb <- merged_df$length_bp / 1000  # Convert bp to kb

# --- CPM ---
cpm <- t(t(counts) / colSums(counts) * 1e6)

# --- RPKM ---
rpkm <- t(t(counts) / colSums(counts) / gene_length_kb * 1e3)

# --- TPM ---
# Step 1: RPK
rpk <- counts / gene_length_kb

# Step 2: Sum RPK per sample
tpm <- t(t(rpk) / colSums(rpk) * 1e6)

# Combine results for inspection
cbind(merged_df[,1:2], CPM_conA = cpm[,1], RPKM_conA = rpkm[,1], TPM_conA = tpm[,1])

print(head(rpkm))
print(head(cpm))
print(head(tpm))

# Add gene IDs as the first column
rpkm_df <- cbind(gene_id = merged_df$gene_id, as.data.frame(rpkm))
cpm_df  <- cbind(gene_id = merged_df$gene_id, as.data.frame(cpm))
tpm_df  <- cbind(gene_id = merged_df$gene_id, as.data.frame(tpm))

# Print
head(rpkm_df)
head(cpm_df)
head(tpm_df)
