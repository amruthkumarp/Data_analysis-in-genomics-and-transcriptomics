library(edgeR)
library(ggplot2)
library(pheatmap)
# Data import 
counts <- read.table("all_samples_counts.csv",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)



# Make metadata
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "treated", "treated")
)



# Make DGEList
group <- factor(coldata$condition)  # conditions (e.g. control, treated)
y <- DGEList(counts=counts, group=group)

# 3. Filter lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,, keep.lib.sizes=FALSE]

# 4. Normalize
y <- calcNormFactors(y)

# 5. Estimate dispersion
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# 6. Fit model and test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)   # compare treated vs control
res <- topTags(lrt, n=Inf)$table
res
res$gene <- rownames(res)

plotMD(lrt, column=1,
       main="MA Plot (edgeR)",
       xlab="Average logCPM", ylab="Log2 Fold Change")
abline(h=c(-1,1), col="blue", lty=2)


# -------------------------------
# 🔥 Volcano plot
# -------------------------------
res$threshold <- as.factor(
  ifelse(res$FDR < 0.05 & abs(res$logFC) >= 1, "Significant", "Not Significant")
)

ggplot(res, aes(x=logFC, y=-log10(FDR), color=threshold)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 FDR") +
  ggtitle("Volcano Plot (edgeR)")

topGenes <- res[res$FDR < 0.05 & abs(res$logFC) >= 1, ]
head(topGenes[order(topGenes$FDR), ], 20)


# -------------------------------
# 🔥 Heatmap of top variable genes
# -------------------------------
# Use logCPM values for visualization
logCPM <- cpm(y, log=TRUE, prior.count=1)

# Pick top 50 most variable genes
topVarGenes <- head(order(apply(logCPM, 1, var), decreasing=TRUE), 50)
mat <- logCPM[topVarGenes, ]

# Z-score transformation (row-wise scaling)
mat <- t(scale(t(mat)))

# Heatmap with clustering
pheatmap(mat,
         annotation_col = coldata,  # sample annotation
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6,
         main = "Heatmap of Top 50 Variable Genes (edgeR)")

