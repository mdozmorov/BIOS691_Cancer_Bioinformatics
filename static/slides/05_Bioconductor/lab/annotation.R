library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(KEGGREST)
library(AnnotationHub)
library(GenomicRanges)
library(ExperimentHub)
library(SummarizedExperiment)
library(dplyr)

# Annotation resources
?org.Hs.eg.db
columns(org.Hs.eg.db)
mapIds(org.Hs.eg.db, keys = c("BRCA1", "BRCA2"), column = "ENSEMBL", keytype = "SYMBOL")

# Extract genomic features from an object
exons(TxDb.Hsapiens.UCSC.hg38.knownGene)
?exonsBy
exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")

# Annotations using BiomaRt
## Discover and then selected mart
head(listMarts(), 7) ## list the marts
head(listDatasets(useMart("ensembl")), 7) ## mart datasets
ensembl <- ## fully specified mart
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")

head(listFilters(ensembl), 7) ## filters
listFilters(ensembl)[grep("entrez", listFilters(ensembl)[, "description"], ignore.case = TRUE), ]

myFilter <- "chromosome_name"
substr(filterOptions(myFilter, ensembl), 1, 50) ## return values
myValues <- c("21", "22")
head(listAttributes(ensembl), 3) ## attributes
myAttributes <- c("ensembl_gene_id", "chromosome_name")

## assemble and query the mart
res <- getBM(
  attributes = myAttributes, filters = myFilter,
  values = myValues, mart = ensembl
)
res

# Get GO annotations for gene names
getBM(attributes = c("hgnc_symbol", "go_id"), filters = "hgnc_symbol", values = c("MYCN"), mart = ensembl)


# AnnotationHub
ah <- AnnotationHub()
ah
unique(ah$dataprovider)
# Subsetting
ah <- subset(ah, species == "Homo sapiens")
ah
query(ah, "grasp2") # see library(grasp2db)
# res <- ah[["AH21414"]] # download actual data
ah_query <- query(ah, c("ucsc", "hg19"))
ah_query
# mcols(ah_query) %>% View()
display(ah)
# Quick look at the specific fields
ah_query$title
ah_query$genome
ah_query$tags
# File types available
ah_query$description %>%
  table() %>%
  sort(., decreasing = TRUE) %>%
  head()
ucsc_data <- query(ah, c("UCSC"))
ucsc_data$title
ucsc_data <- query(ah, c("UCSC", "GWAS Catalog"))
ucsc_data
ucsc_data$genome
ucsc_data <- query(ah, c("UCSC", "GWAS Catalog", "hg19"))
ucsc_data

gwas_granges <- ah_query[["AH5029"]]
gwas_granges
summary(width(gwas_granges))
summary(width(reduce(gwas_granges)))
#
ucsc_selected <- mcols(ah_query)
ucsc_selected
ucsc_selected <- ucsc_selected[grepl("GRanges object from UCSC track", ucsc_selected$description) & !grepl("Chain", ucsc_selected$description), ]
ucsc_selected
chromosome_bands <- ah_query[["AH5012"]]
chromosome_bands
metadata(chromosome_bands)
chromosome_bands@elementMetadata
query(AnnotationHub(), c("release-88", "homo"))

library(ExperimentHub)
eh <- ExperimentHub()
eh
head(unique(eh$dataprovider))
head(unique(eh$species))
tcga_data <- query(ExperimentHub(), c("TCGA"))
tcga_data$title[grepl("BRCA", tcga_data$title)]
tcga_data <- query(ExperimentHub(), c("BRCA_RNASeqGene-20160128"))
tcga_data
brca_data <- tcga_data[["EH597"]]
brca_data
assay(brca_data)[1:5, 1:5]
