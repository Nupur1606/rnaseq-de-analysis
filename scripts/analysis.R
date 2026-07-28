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