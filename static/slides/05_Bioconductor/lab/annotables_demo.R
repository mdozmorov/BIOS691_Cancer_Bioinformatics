library(annotables)
# Remove non-canonical chromosome names
grch37 <- grch37[ !(grepl("_", grch37$chr) | grepl("GL", grch37$chr)), ]
# Replace "MT" by "M"
grch37$chr <- gsub("MT", "M", grch37$chr)
# Append "chr" prefix
grch37$chr <- paste("chr", grch37$chr, sep="")
# Replace missing gene names and EntrezIDs by "?"
grch37$entrez[ is.na(grch37$entrez) ] <- "?"
grch37$symbol[ is.na(grch37$symbol) ] <- "?"
# Replace strand
grch37$strand[ grch37$strand == -1] <- "-"
grch37$strand[ grch37$strand ==  1] <- "+"

# Exploring how many transcript types we have
# Description of different biotypes
# http://www.ensembl.org/Help/Faq?id=468
types <- grch37$biotype %>% table() %>% data.frame() # Get counts of unique ones
types <- types[order(types$Freq, decreasing = T), ] # Sort from highest to lowest counts

# Extract coordinates
for (t in types$.) {
  print(t)
  mtx <- grch37[ grch37$biotype == t, , drop = FALSE]
  mtx <- data.frame(mtx$chr, mtx$start, mtx$end, name=paste(mtx$entrez, mtx$symbol, sep = "|"), mtx$strand)
  write.table(mtx, paste(t, ".bed", sep=""), sep="\t", quote = F, col.names = F, row.names = F)
}
