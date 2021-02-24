# Retrieve the coding sequences grouped by transcript for the gene of interest and verify that each coding sequence is a multiple of 3.

library(EnsDb.Hsapiens.v86) # hg19/GRCh37 coordinates
edb <- EnsDb.Hsapiens.v86

# Use the `cdsBy()` function to retrieve the genomic coordinates of all coding
# sequences for the gene 'BRCA1' from the [EnsDb.Hsapiens.v86][] package. To
# retrieve only data for the specified gene, submit either a `GenenameFilter`
# or a filter formula/expression to the function's `filter` parameter. This
# avoids to extract the coding region for all genes, which takes a long time.

brca1cds <- cdsBy(edb, by = "tx", filter = ~ genename == "BRCA1")

class(brca1cds)
length(brca1cds)
brca1cds[[1]]                           # exons in cds
cdswidth <- width(brca1cds)             # width of each exon
all((sum(cdswidth) %% 3) == 0)          # sum within cds, modulus 3

# The CDS for some transcripts is not of the expected length, how comes? Get the
# transcript ID of the first transcript that does have a CDS of the wrong size and
# look this transcript up in the Ensembl genome browser
# (http://www.ensembl.org).

tx_cds_fail <- names(brca1cds)[(sum(cdswidth) %% 3) != 0]
length(tx_cds_fail)
tx_cds_fail[1]

# In the description of the transcript it says *CDS 5' incomplete*. Thus, in
# addition to known protein coding transcripts, Ensembl provides also annotations
# for transcripts known to be targeted for nonsense mediated mRNA decay or that
# have incomplete CDS. Such transcripts would however not be listed in e.g.
# the [TxDb.Hsapiens.UCSC.hg19.knownGene][] package.

# Next we visualize the BRCA1 transcripts using [Gviz][] (this package has an
# excellent vignette, `vignette("Gviz")`)

library(Gviz)
## Use the function from the ensembldb package to extract the data in the
## format suitable for Gviz
grt <- getGeneRegionTrackForGviz(edb, filter = ~genename == "BRCA1")
plotTracks(list(GenomeAxisTrack(), GeneRegionTrack(grt)))

# Extract the coding sequences of each transcript. `EnsDb` databases provide
# annotations from Ensembl and use hence Ensembl style chromosome names (such as
# "Y") while the `BSgenome` package is based on UCSC annotations that use a naming
# style that prepends a "chr" to each chromosome name (e.g. "chrY"). Change thus
# the `seqlevelsStyle` from the default UCSC chromosome naming to Ensembl naming
# style.

library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

## Change the seqlevelsStyle from UCSC to Ensembl
seqlevelsStyle(genome) <- "Ensembl"
tx_seq <- extractTranscriptSeqs(genome, brca1cds)
tx_seq

# We can also inspect the CDS sequence for the transcripts with incomplete
# CDS. Many of them do not start with a start codon hence indicating that the CDS
# is incomplete on their 5' end.

tx_seq[tx_cds_fail]

# Intron coordinates can be identified by first calculating the range of
# the genome (from the start of the first exon to the end of the last
#             exon) covered by each transcript, and then taking the (algebraic) set
# difference between this and the genomic coordinates covered by each
# exon

introns <- psetdiff(unlist(range(brca1cds)), brca1cds)

# Retrieve the intronic sequences with `getSeq()` (these are *not*
#                                                    assembled, the way that `extractTranscriptSeqs()` assembles exon
#                                                  sequences into mature transcripts); note that introns start and end
# with the appropriate acceptor and donor site sequences.
# Unfortunately, UCSC and Ensembl do also use different names for the genome
# assembly. Change the genome name for the `introns` object to matche the one from
# the `genome` object.

unique(genome(genome))
genome(introns)

## Change the genome name on introns to macht the one from the
## BSgenome package
genome(introns) <- c(`17` = unique(genome(genome)))

seq <- getSeq(genome, introns)
names(seq)
seq[["ENST00000352993"]]                     # 20 introns
```

# Here we use [rtracklayer][] to retrieve estrogen receptor binding
# sites identified across cell lines in the ENCODE project. We focus on
# binding sites in the vicinity of a particularly interesting region of
# interest.
# 
# 1. Define our region of interest by creating a `GRanges` instance with
# appropriate genomic coordinates. Our region corresponds to 10Mb up-
#   and down-stream of a particular gene.
# 2. Create a session for the UCSC genome browser
# 3. Query the UCSC genome browser for ENCODE estrogen receptor
# ERalpha<sub>a</sub> transcription marks; identifying the
# appropriate track, table, and transcription factor requires
# biological knowledge and detective work.
# 4. Visualize the location of the binding sites and their scores;
# annotate the mid-point of the region of interest.

# Define the region of interest
library(GenomicRanges)
roi <- GRanges("chr10", IRanges(92106877, 112106876, names="ENSG00000099194"))

# Create a session
library(rtracklayer) 
session <- browserSession()
genome(session) <- "hg19"
tracks <- trackNames(session)
tracks[grep("tf", tracks, ignore.case = TRUE)]

# Query the UCSC for a particular track, table, and transcription
# factor, in our region of interest

trackName <- "wgEncodeTfBindingSuper"
tableName <- "wgEncodeTfBindingSuper"
trFactor <- "ERalpha_a"

query <- ucscTableQuery(session, range = roi)
# track(query)
tables <- tableNames(query)
tables[grep("tfbsclustered", tables, ignore.case = TRUE)]
tableName(query) <- "wgEncodeRegTfbsClusteredV2"

names(query) <- trFactor
ucscTable <- getTable(query)

# Visualize the result
plot(score ~ chromStart, ucscTable, pch="+")
abline(v=start(roi) + (end(roi) - start(roi) + 1) / 2, col="blue")
