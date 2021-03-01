rm(list = ls())
library(biomaRt)
# Some biomaRt hosts host='caprica.caltech.edu' host='www.sanger.ac.uk' host='biomart.informatics.jax.org'
# host='biomart.intogen.org' host='ensembl.gramene.org'

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Archived biomart version for mm9 mart=useMart(host='may2009.archive.ensembl.org',
# biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl') or the latest biomart version for mm mart <-
# useMart('ensembl', dataset='mmusculus_gene_ensembl') # rnorvegicus or the latest for hs mart <-
# useMart('ensembl', dataset='hsapiens_gene_ensembl', host='www.biomart.org')

# Information - lists of filters and attributes
db_list <- listDatasets(mart) # List available datasets
filters <- listFilters(mart) # Source IDs to convert
head(filters)
filters[grep("illumina", filters[, "name"]), ] # List specific ones
attr <- listAttributes(mart) # Destination IDs
head(attr)
attr[grep("symbol", attr[, 1]), ] # List gene symbol-specific
attr[grep("go_", attr[, 1]), ] # List GO-specific

IDs <- readLines("clipboard") # Read a list of Ensembl IDs
IDs <- c("ENSG00000205420", "ENSG00000206075")
# Convert from Ensembl IDs to Gene Name and Description
genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), filters = "ensembl_gene_id", values = IDs,
  mart = mart
) # , uniqueRows=T)
# Convert from UCSC IDs to Gene Name and Description
IDs <- "uc001aak.4"
genes <- getBM(attributes = c("ucsc", "hgnc_symbol", "description"), filters = "ucsc", values = IDs, mart = mart) # , uniqueRows=T)
# Convert Gene names IDs (hgnc_symbol, mgi_symbol)
IDs <- c("AKT1", "IFNG", "STAT1")
genes <- getBM(
  attributes = c("wikigene_name", "hgnc_symbol", "description", "ensembl_gene_id"), filters = "hgnc_symbol",
  values = IDs, mart = mart
) # , uniqueRows=T)
# Convert entrezgene
genes <- getBM(attributes = c("hgnc_symbol", "description"), filters = "entrezgene", values = "1", mart = mart, uniqueRows = T)
# Convert RefSEq IDs
IDs <- c("NM_001006946")
genes <- getBM(
  attributes = c("refseq_mrna", "hgnc_symbol", "description"), filters = "refseq_mrna", values = IDs,
  mart = mart
) # , uniqueRows=T)
# Get GO annotations for gene names
IDs <- c("AKT1", "IFNG", "STAT1")
getBM(attributes = c("hgnc_symbol", "go_id"), filters = "hgnc_symbol", values = IDs, mart = mart)

write.table(genes, "clipboard-128", sep = "\t")
write.table(genes, "F:/111.txt", sep = "\t")


library("AnnotationDbi")
library("org.Hs.eg.db")
res$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = " rst")
res$entrez <- mapIds(org.Hs.eg.db, keys = row.names(res), column = "ENTREZID", keytype = "ENSEMBL", multiVals = " rst")
