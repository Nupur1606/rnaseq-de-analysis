# ================================================
# RNA-seq Differential Expression Analysis
# Dataset: Airway (dexamethasone treated vs untreated)
# ================================================

# Load libraries
library(SummarizedExperiment)
library(DESeq2)
library(airway)

# Custom counts function (workaround for S4 method issue)
counts <- function(x) x@assays@data$counts

# Load dataset
data(airway)

# Explore the data
dim(airway)           # how many genes and samples
colData(airway)       # sample metadata
head(counts(airway))  # first 6 rows of count matrix

# ================================================
# Step 2: Filter low expression genes
# ================================================

# Keep only genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts(airway) >= 10) >= 3
airway_filtered <- airway[keep, ]

# How many genes remain after filtering?
dim(airway_filtered)

# ================================================
# Step 3: Run DESeq2
# ================================================

# Tell DESeq2 which column defines our comparison
dds <- DESeqDataSet(airway_filtered, design = ~ cell + dex)

# Run the full DESeq2 pipeline
dds <- DESeq(dds)

# ================================================
# Step 4: Extract results
# ================================================

# Extract results for treated vs untreated comparison
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Summary of results
summary(res)

# Look at first few rows
head(res)

# ================================================
# Step 5: Save results
# ================================================

# Convert to dataframe and save
res_df <- as.data.frame(res)
write.csv(res_df, 
          file = "../results/deseq2_results.csv",
          row.names = TRUE)

cat("Results saved! Total significant genes (padj < 0.05):", 
    sum(res_df$padj < 0.05, na.rm = TRUE), "\n")

    # ================================================
# Step 6: Volcano Plot
# ================================================

library(ggplot2)

# Prepare data for plotting
res_df$significant <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant", "Not Significant"
)

res_df$significant[is.na(res_df$significant)] <- "Not Significant"

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, 
                                    y = -log10(padj),
                                    color = significant)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c("Not Significant" = "grey60",
                                 "Significant" = "#E63946")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", 
             color = "black", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", linewidth = 0.3) +
  labs(title = "Volcano Plot: Dexamethasone Treated vs Untreated",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value",
       color = "Significance") +
  theme_minimal()
  # Save plot
ggsave("../results/volcano_plot.png", 
       plot = volcano_plot,
       width = 8, height = 6, dpi = 300)

cat("Volcano plot saved!\n")

# ================================================
# Step 7: Heatmap of top 50 DE genes
# ================================================

library(pheatmap)
library(RColorBrewer)

# Get top 50 significant genes by adjusted p-value
top_genes <- rownames(res_df[order(res_df$padj), ])[1:50]

# Extract normalized counts for these genes
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, ]

# Scale each gene (row) so colors show relative expression
mat_scaled <- t(scale(t(mat)))

# Create annotation for columns (samples)
annotation <- data.frame(
  Treatment = colData(dds)$dex,
  row.names = colnames(mat)
)

# Plot heatmap
pheatmap(mat_scaled,
         annotation_col = annotation,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes",
         filename = "../results/heatmap.png",
         width = 8,
         height = 10)

cat("Heatmap saved!\n")