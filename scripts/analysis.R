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