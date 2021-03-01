# Summarized Experiments
dataDir <- "/Users/mdozmorov/Documents/Work/Teaching/BIOS691_Cancer_Bioinformatics/static/slides/05_Bioconductor/lab/airway_data"

# three components -- underlying 'matrix', 'row' annotations (genomic features), 'column' annotations (sample descriptions)
counts <- read.csv(file.path(dataDir, "airway_counts.csv"), row.names=1)
counts <- as.matrix(counts)
counts[1:5, 1:5]
dim(counts)

colData <- read.csv(file.path(dataDir, "airway_colData.csv"), row.names=1)
colData[, 1:4]

rowRanges <- readRDS(file.path(dataDir, "airway_rowRanges.rds"))
rowRanges

# Plotting counts in treated vs. untreated groups - easy to make mistakes in ordering
cidx <- colData$dex == "trt"
plot(
  rowMeans(1 + counts[, cidx]) ~ rowMeans(1 + counts[, !cidx]),
  log="xy"
)

# Making one object that will have all parts aligned with each other, and accessible
library(SummarizedExperiment)

se <- SummarizedExperiment(counts, 
                           rowRanges = rowRanges, 
                           colData = colData)
se
se$dex
colData(se)[, "dex"]
cidx <- se$dex == "trt"
plot(
  rowMeans(1 + assay(se)[, cidx]) ~ rowMeans(1 + assay(se)[, !cidx]),
  log="xy"
)

# Much more robust to 'bookkeeping' errors
# matrix-like interface: `dim()`, two-dimensional `[`, ...
# accessors: `assay()`, `rowData()` / `rowRanges()`, `colData()`, ...

methods(class="SummarizedExperiment")
?"SummarizedExperiment"
browseVignettes(package="SummarizedExperiment")
