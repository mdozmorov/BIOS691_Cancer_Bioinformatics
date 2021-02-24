- `01_Bioconductor` - introduction to Bioconductor, main genomic packages
- `Annotation.Rmd` - gene IDs, biomaRt

- `biostrings.pdf` - Biostrings lecture. http://www.haowulab.org/teaching/bioc/biostrings.pdf
- `lab/biostrings.R` - R code

- `genomicRanges.pdf` - GenomicRanges and GenomicFeatures, IRanges/GRanges, TxDb, biomaRt. `lab/genomicRanges.R` - R code https://github.com/haowulab/haowulab.github.io/blob/master/teaching/bioc/genomicRanges.pdf

- `lab/annotation.R` - OrgDb, TxDb, biomaRt, AnnotationHub, ExperimentHub

- `lab/shortread.R` - basics of ShortRead package, read in FASTQ files, getting QC. Uses `fastq` data

- `lab/SummarizedExperiment.R` - SummarizedExperiment using 'airway' data.


## `lab` - in-class lab material

- `ranges.Rmd` - IRanges/GRanges, ensembldb, gene length distribution, overlap on GRanges
- `GRL.Rmd` - lists of GRanges, EnsDb.Hsapiens.v86, Gviz visualization
- `SE.Rmd` - SummarizedExperiment demo by Michael Love, manual construction

- `ensembl_gviz_rtracklayer.R` - Gviz visualization
- `Install_Bioconductor.Rmd` - Installation instructions
- `eSet.Rmd` - ExpressionSet, annotations, SummarizedExperiment
- `Annotation.Rmd` - working with annotation packages
- `GEO.Rmd` - working with GEO
- `biomaRt.R` - Biomart tutorial
- `annotables_demo.R` - annotables demo, extracting gene classes


## References

- Gentleman, Robert C., Vincent J. Carey, Douglas M. Bates, Ben Bolstad, Marcel Dettling, Sandrine Dudoit, Byron Ellis, et al. “[Bioconductor: Open Software Development for Computational Biology and Bioinformatics](https://doi.org/10.1186/gb-2004-5-10-r80).” Genome Biology 5, no. 10 (2004): R80. 

- Huber, Wolfgang, Vincent J. Carey, Robert Gentleman, Simon Anders, Marc Carlson, Benilton S. Carvalho, Hector Corrada Bravo, et al. “[Orchestrating High-Throughput Genomic Analysis with Bioconductor](https://doi.org/10.1038/nmeth.3252).” Nature Methods, (February 2015) - Bioconductor overview

- Ramos, Marcel, Lucas Schiffer, Angela Re, Rimsha Azhar, Azfar Basunia, Carmen Rodriguez Cabrera, Tiffany Chan, et al. “[Software For The Integration Of Multi-Omics Experiments In Bioconductor](http://dx.doi.org/10.1158/0008-5472.can-17-0344),” June 1, 2017. - MultiAssayExperiment, Bioconductor package for management of multi-assay data

- Lawrence, Michael, Wolfgang Huber, Hervé Pagès, Patrick Aboyoun, Marc Carlson, Robert Gentleman, Martin T. Morgan, and Vincent J. Carey. “[Software for Computing and Annotating Genomic Ranges](https://doi.org/10.1371/journal.pcbi.1003118).” PLoS Computational Biology, (2013) - IRanges, GenomicRanges, GenomicFeatures packages.

### Bioconductor resources

- "[Bioconductor for Genomics Data Science](https://kasperdanielhansen.github.io/genbioconductor/)" course by Kasper Hansen. IRanges, GRanges, seqinfo, AnnotationHub, Biostrings, BSgenome, rtracklayer, and more. Includes videos, code examples and lecture material. [GitHub](https://github.com/kasperdanielhansen/genbioconductor)
 
- [Introduction to Computational Biology](https://biodatascience.github.io/compbio/) course by Mike Love. [GitHub](https://github.com/biodatascience/compbio)
 
- [Multi-omic Integration and Analysis of cBioPortal and TCGA data with MultiAssayExperiment](https://github.com/waldronlab/MultiAssayWorkshop) - workshop by Levi Waldron's group
 

# `misc`

- `TutorialDurinck.ppt` - detailed presentation of Bioconductor, analysis, [http://www.nettab.org/2003/docs/TutorialDurinck.ppt](http://www.nettab.org/2003/docs/TutorialDurinck.ppt)

# ToDo

- Kishore, Kamal, Stefano de Pretis, Ryan Lister, Marco J. Morelli, Valerio Bianchi, Bruno Amati, Joseph R. Ecker, and Mattia Pelizzola. “MethylPipe and CompEpiTools: A Suite of R Packages for the Integrative Analysis of Epigenomics Data.” BMC Bioinformatics 16 (September 29, 2015): 313. https://doi.org/10.1186/s12859-015-0742-6. - methylPipe + compEpiTools. Integrative analysis of methylation and other omics data. Detailed description of each function. Table 1 lists comparison with other tools. https://bioconductor.org/packages/release/bioc/html/methylPipe.html, https://bioconductor.org/packages/release/bioc/html/compEpiTools.html, https://bioconductor.org/packages/release/data/experiment/html/ListerEtAlBSseq.html

- Gu, Zuguang, Roland Eils, Matthias Schlesner, and Naveed Ishaque. “EnrichedHeatmap: An R/Bioconductor Package for Comprehensive Visualization of Genomic Signal Associations.” BMC Genomics 19, no. 1 (April 2018): 234. https://doi.org/10.1186/s12864-018-4625-x. - EnrichedHeatmap - genomics data visualization, intro on other packages (deeptools, ngs.plot, genomation). https://bioconductor.org/packages/release/bioc/html/EnrichedHeatmap.html

- Argelaguet, Ricard, Britta Velten, Damien Arnol, Sascha Dietrich, Thorsten Zenz, John C Marioni, Florian Buettner, Wolfgang Huber, and Oliver Stegle. “Multi‐Omics Factor Analysis—a Framework for Unsupervised Integration of Multi‐omics Data Sets.” Molecular Systems Biology 14, no. 6 (June 2018): e8124. https://doi.org/10.15252/msb.20178124. - MOFA - Multi-Omics Factor Analysis, infers factors that capture biological and technical variability. Separates variability shared and specific across datasets. Used for identification of sample subgroups, imputation, detection of outliers. Refs to multi-omics single-cell technologies. Group Factor Analysis-based (Fig 1), good methods (supplemental) and interpretation description.  Compared with iCluster and canonical GFA. MOFAtools R package interfaced with a Python module. Tutorials at https://github.com/bioFAM/MOFA

- List of software packages for multi-omics analysis. https://github.com/mikelove/awesome-multi-omics

- Huang, Sijia, Kumardeep Chaudhary, and Lana X. Garmire. “More Is Better: Recent Progress in Multi-Omics Data Integration Methods.” Frontiers in Genetics 8 (2017): 84. https://doi.org/10.3389/fgene.2017.00084. - Review of multi-omics integration methods and individual tools. Unsupervised integration (Matrix Factorization, correlation-based, Bayesian methods, network-based methods, multiple kernel learning). Supervised integration (network-based methods, multiple kernel learning). Semi-supevised methods. Biological insights, survival predictions.
