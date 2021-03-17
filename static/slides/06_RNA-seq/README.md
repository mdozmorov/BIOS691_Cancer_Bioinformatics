

## `lab` - in-class lab material

- `rnaseqGene_CSAMA2017.Rmd` - "RNA-seq workflow: gene-level exploratory analysis and differential expression," https://www.bioconductor.org/help/course-materials/2017/CSAMA/labs/2-tuesday/lab-03-rnaseq/rnaseqGene_CSAMA2017.html

- `rnaseqGene` - RNA-seq workflow: gene-level exploratory analysis and differential expression. https://github.com/mikelove/rnaseqGene

## `read` - reading assignment

### Overview

- Lowe, Rohan, Neil Shirley, Mark Bleackley, Stephen Dolan, and Thomas Shafee. “[Transcriptomics Technologies](https://doi.org/10.1371/journal.pcbi.1005457).” PLoS Computational Biology, (May 2017)

- Conesa, Ana, Pedro Madrigal, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, et al. “[A Survey of Best Practices for RNA-Seq Data Analysis](https://doi.org/10.1186/s13059-016-0881-8).” Genome Biology, (December 2016) - RNA-seq analysis roadmap, QC. Differential detection. TPM. Tools for alternative splicing detection and visualization. small RNA analysis. single cell. Integrative analysis, with methylation.

- Garber, Manuel, Manfred G. Grabherr, Mitchell Guttman, and Cole Trapnell. “[Computational Methods for Transcriptome Annotation and Quantification Using RNA-Seq](https://doi.org/10.1038/nmeth.1613).” Nature Methods, (June 2011) - RNA-seq alignment and quantification. Table of tools. Transcriptome reconstruction. Alternative splicing

- Wang, Zhong, Mark Gerstein, and Michael Snyder. “[RNA-Seq: A Revolutionary Tool for Transcriptomics](https://doi.org/10.1038/nrg2484).” Nature Reviews. Genetics, (January 2009) - RNA-seq review. 

- Williams, Alexander G., Sean Thomas, Stacia K. Wyman, and Alisha K. Holloway. “[RNA-Seq Data: Challenges in and Recommendations for Experimental Design and Analysis: RNA-Seq Data: Experimental Design and Analysis]((https://doi.org/10.1002/0471142905.hg1113s83)).” In Current Protocols in Human Genetics, John Wiley & Sons, Inc., 2014. - RNA-seq basics, tools, simulations

- Altman, Naomi, and Martin Krzywinski. “[Points of Significance: Sources of Variation](https://doi.org/10.1038/nmeth.3224).” Nature Methods, (December 30, 2014)

- Martin, Jeffrey A., and Zhong Wang. “[Next-Generation Transcriptome Assembly](https://doi.org/10.1038/nrg3068).” Nature Reviews. Genetics, (September 7, 2011) - Transcriptome assembly. Sequencing technologies overview. Reference-based and de novo assembly, combined approach idea. Splice graph, De Bruijn graph.

- Marioni, John C., Christopher E. Mason, Shrikant M. Mane, Matthew Stephens, and Yoav Gilad. “[RNA-Seq: An Assessment of Technical Reproducibility and Comparison with Gene Expression Arrays](https://doi.org/10.1101/gr.079558.108).” Genome Research, (September 2008) - Illumina sequencing - microarray comparison. Good agreement. Assessing lane effect with hypergeometric distribution. Likelihood ratio test for differential expression. Chi-squared goodness-of-fit test.

- Peixoto, Lucia, Davide Risso, Shane G. Poplawski, Mathieu E. Wimmer, Terence P. Speed, Marcelo A. Wood, and Ted Abel. “[How Data Analysis Affects Power, Reproducibility and Biological Insight of RNA-Seq Studies in Complex Datasets](https://doi.org/10.1093/nar/gkv736).” Nucleic Acids Research, (September 18, 2015) - The importance of RNA-seq normalization and batch effect removal. RUVseq increases power, but the choice of the number of latent variables is important. [Tutorial: Steps in RNA-seq data processing, normalization, exploratory data analysis](https://github.com/drisso/peixoto2015_tutorial)

- Wang, Eric T., Rickard Sandberg, Shujun Luo, Irina Khrebtukova, Lu Zhang, Christine Mayr, Stephen F. Kingsmore, Gary P. Schroth, and Christopher B. Burge. “[Alternative Isoform Regulation in Human Tissue Transcriptomes](https://doi.org/10.1038/nature07509).” Nature, (November 27, 2008) - Alternative splicing comparison between tissues. \~94% of genes are alternatively transcribed. Variation in alternative splicing is much more between tissues than between individuals.

- Park, Eddie, Zhicheng Pan, Zijun Zhang, Lan Lin, and Yi Xing. “[The Expanding Landscape of Alternative Splicing Variation in Human Populations](https://doi.org/10.1016/j.ajhg.2017.11.002).” The American Journal of Human Genetics, (January 2018) - Alternative splicing, detailed overview

### Statistics

- Pachter, Lior. “[Models for Transcript Quantification from RNA-Seq](https://arxiv.org/abs/1104.3889).” ArXiv, 2011. - RNA-seq quantification statistics, expectation-maximization algorithm

- Robinson, Mark D., and Alicia Oshlack. “[A Scaling Normalization Method for Differential Expression Analysis of RNA-Seq Data](https://doi.org/10.1186/gb-2010-11-3-r25).” Genome Biology, (2010) - TMM normalization method. Problems with library scaling normalization. Well-written intuitive motivating example. MA plot, trimming outliers, weighted (inverse of the variance) M average after discarding 30% of M outliers and lowest 5% of A values.

- "[How not to perform a differential expression analysis (or science)](https://liorpachter.wordpress.com/2017/08/02/how-not-to-perform-a-differential-expression-analysis-or-science/)" blog post by Lior Pachter, about Salmon-kallisto similarities and differences, general references

- Robinson, Mark D., and Gordon K. Smyth. “[Small-Sample Estimation of Negative Binomial Dispersion, with Applications to SAGE Data](https://doi.org/10.1093/biostatistics/kxm030).” Biostatistics (Oxford, England), (April 2008) - Negative Binomial distribution instead of Poisson. Previous models: binomial, Poisson.

- Lun, Aaron T. L., and Gordon K. Smyth. “[No Counts, No Variance: Allowing for Loss of Degrees of Freedom When Assessing Biological Variability from RNA-Seq Data](https://doi.org/10.1515/sagmb-2017-0010).” Statistical Applications in Genetics and Molecular Biology, (April 25, 2017) - Negative impact of genes with zero counts on GLM framework for RNA-seq differential expression analysis. Overdispersion, GLM, quasi-likelihood F-test, adjusting degrees of freedom for zero-count genes.

- Law, Charity W, Yunshun Chen, Wei Shi, and Gordon K Smyth. “[Voom: Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read Counts](https://doi.org/10.1186/gb-2014-15-2-r29).” Genome Biology, (2014) - voom paper

- Love, Michael I, Wolfgang Huber, and Simon Anders. “[Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8).” Genome Biology, (December 2014) - DESeq2 paper. Problems with fold-change ranking of genes - proposed solution using shrinkage of FCs. Generalized linear model using Negative Binomial distribution. Borrowing information - genes of similar avarage expression have similar dispersion. rlog-transformation. [The original DESeq publication](https://www.ncbi.nlm.nih.gov/pubmed/20979621)

- Witten, Daniela M. “[Classification and Clustering of Sequencing Data Using a Poisson Model](https://doi.org/10.1214/11-AOAS493).” The Annals of Applied Statistics, (December 2011) - RNA-seq modeling with Poisson distribution. samples X genes matrix. Derivation of Poisson, negative binomial, using Poisson for linear discriminant analysis and clustering (Poisson dissimilarity).

- Patro, Rob, Geet Duggal, Michael I Love, Rafael A Irizarry, and Carl Kingsford. “[Salmon Provides Fast and Bias-Aware Quantification of Transcript Expression](https://doi.org/10.1038/nmeth.4197).” Nature Methods, (March 6, 2017) - Salmon paper. Pseudo-alignment, or using precomputed alignment to tramscriptome. Dual-phase statistical inference procedure and sample-specific bias models that account for sequence-specific, fragment, GC content, and positional biases. Comparison with kallisto and sailfish. Tests on simulated (Polyester, RSEM-sim) and real (GEUVADIS, SEQC) data. Detailed Methods description. [COMBINE-lab/Salmon GitHub repo](https://github.com/COMBINE-lab/Salmon)

- Jiang, Hui, and Wing Hung Wong. “[Statistical Inferences for Isoform Expression in RNA-Seq](https://doi.org/10.1093/bioinformatics/btp113).” Bioinformatics, (April 15, 2009) - Alternative splicing statistics - Poisson modeling. Problem - most reads are shared by more than one isoform. How to quantify isoform expression from exon counts. Detailed statistical derivations

- Young, Matthew D., Matthew J. Wakefield, Gordon K. Smyth, and Alicia Oshlack. “[Gene Ontology Analysis for RNA-Seq: Accounting for Selection Bias](https://doi.org/10.1186/gb-2010-11-2-r14).” Genome Biology, (2010) - Gene set enrichment analysis accounting for length and expression of transcripts. Instead of random sampling, use of the Wallenius non-central hypergeometric distribution to account for biased sampling. [GOseq R package](https://bioconductor.org/packages/goseq/)

- Law, Charity W., Kathleen Zeglinski, Xueyi Dong, Monther Alhamdoosh, Gordon K. Smyth, and Matthew E. Ritchie. “[A Guide to Creating Design Matrices for Gene Expression Experiments.](https://doi.org/10.12688/f1000research.27893.1)” F1000Research (December 10, 2020) - Design matrices for various experimental designs. Means model or mean-reference model.

- Soneson, C, F Marini, F Geier, MI Love, and MB Stadler. “[ExploreModelMatrix: Interactive Exploration for Improved Understanding of Design Matrices and Linear Models in R](https://doi.org/10.12688/f1000research.24187.2)” F1000Research, (June 4, 2020). 

### Workflows and tools

- He, Wen, Shanrong Zhao, Chi Zhang, Michael S. Vincent, and Baohong Zhang. “[QuickRNASeq: Guide for Pipeline Implementation and for Interactive Results Visualization](https://doi.org/10.1007/978-1-4939-7710-9_4).” Springer New York, 2018. - Practical RNA-seq tutorial based on QuickRNASeq publication.

- Love, Michael I., Simon Anders, Vladislav Kim, and Wolfgang Huber. “[RNA-Seq Workflow: Gene-Level Exploratory Analysis and Differential Expression](https://doi.org/10.12688/f1000research.7035.2).” F1000Research (November 17, 2016) - RNA-seq workflow. From count import, including tximport, through EDA, DESeq2, batch removal, time course analysis, visualization.

- Griffith, Malachi, Jason R. Walker, Nicholas C. Spies, Benjamin J. Ainscough, and Obi L. Griffith. “[Informatics for RNA Sequencing: A Web Resource for Analysis on the Cloud](https://doi.org/10.1371/journal.pcbi.1004393).” PLoS Computational Biology, (August 2015) - RNA-seq technology and analysis introduction. Very interesting are supplementary tables. [The full tutorials](https://rnabio.org/)

- Sahraeian, Sayed Mohammad Ebrahim, Marghoob Mohiyuddin, Robert Sebra, Hagen Tilgner, Pegah T. Afshar, Kin Fai Au, Narges Bani Asadi, et al. “[Gaining Comprehensive Biological Insight into the Transcriptome by Performing a Broad-Spectrum RNA-Seq Analysis](https://doi.org/10.1038/s41467-017-00050-4).” Nature Communications, (December 2017) - RNAcocktail - RNA-seq tools benchmarking. All aspects of RNA-seq analysis, structured, Fig 1. Recommended tools Fig 8. [RNACocktail - A comprehensive framework for accurate and efficient RNA-Seq analysis](https://bioinform.github.io/rnacocktail/), [RNA-seq blog: Unleash the power within RNA-seq](http://www.rna-seqblog.com/unleash-the-power-within-rna-seq/)





- Law, Charity W., Monther Alhamdoosh, Shian Su, Gordon K. Smyth, and Matthew E. Ritchie. “[RNA-Seq Analysis Is Easy as 1-2-3 with Limma, Glimma and EdgeR](https://doi.org/10.12688/f1000research.9005.2).” F1000Research (2016) - [Latest Rsubread-limma plus pipeline](https://f1000research.com/articles/5-1408/v3). [The complete R code for RNA-seq analysis tutorial](https://bioconductor.org/packages/RNAseq123/)

- Pertea, Mihaela, Daehwan Kim, Geo M. Pertea, Jeffrey T. Leek, and Steven L. Salzberg. “[Transcript-Level Expression Analysis of RNA-Seq Experiments with HISAT, StringTie and Ballgown](https://doi.org/10.1038/nprot.2016.095).” Nature Protocols, (September 2016) - New Tuxedo suite. Protocol.

- "[RNA-Seq Methods and Algorithms](https://youtu.be/96yBPM8lEt8)" 7m video by Harold Pimentel, pseudoalignment, kallisto, sleuth, practical

- Zhao, Qi, Yubin Xie, Peng Nie, Rucheng Diao, Licheng Sun, Zhixiang Zuo, and Jian Ren. “[IDEA: A Web Server for Interactive Differential Expression Analysis with R Packages](https://doi.org/10.1101/360461),” July 3, 2018. - Differential expression analysis from a matrix of FPKMs and a design matrix. Several methods to detect DEGs (DESeq2, edgeR, NOISeq, PoissonSeq, SAMseq), plots (MA, volcano, heatmap). [http://renlab.org:3838/IDEA/](http://renlab.org:3838/IDEA/)

- TPMCalculator quantifies mRNA abundance directly from the alignments by parsing BAM files. The input parameters are the same GTF files used to generate the alignments, and one or multiple input BAM file(s) containing either single-end or paired-end sequencing reads. The TPMCalculator output is comprised of four files per sample reporting the TPM values and raw read counts for genes, transcripts, exons and introns respectively. [https://github.com/ncbi/TPMCalculator](https://github.com/ncbi/TPMCalculator)

### RNA-seq resources

- [Tools for RNA-seq data analysis](http://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=info:doi/10.1371/journal.pcbi.1004393.s004)
-[ RNAseq analysis notes by Tommy Tang](https://github.com/crazyhottommy/RNA-seq-analysis)
- [University of Oregon's RNA-seqlopedia](http://rnaseq.uoregon.edu/), a comprehensive guide to RNA-seq starting with experimental design, going through library prep, sequencing, and data analysis.
- RNA-seq blog, [http://www.rna-seqblog.com/](http://www.rna-seqblog.com/), Several blog posts per week on new methods and tools for RNA-seq analysis
- [Informatics for RNA-seq, by Griffith lab](https://github.com/griffithlab/rnaseq_tutorial)

### Practicals

- [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://bioconductor.org/packages/rnaseqGene/) R package

- [RNA-seq analysis exercise using Galaxy](https://usegalaxy.org/u/jeremy/p/galaxy-rna-seq-analysis-exercise), an example analysis using the Tophat+Cufflinks workflow.

- ["bcbio-nextgen" - Validated, scalable, community developed variant calling, RNA-seq and small RNA analysis](https://github.com/chapmanb/bcbio-nextgen)

- [Getting started with HISAT, StringTie, and Ballgown](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/)

- ["enrichOmics" - Functional enrichment analysis of high-throughput omics data](https://github.com/waldronlab/enrichOmics). From basic ExpressionSet differential and functional enrichment analysis to genomic region enrichment analysis and MultiAssayExperiment demo

- [ENCODE RNA-, ChIP-, DNAse-, ATAC- and other seq pipelines](https://github.com/ENCODE-DCC/)



## `misc` - misc presentations and materials

- Video (12 min) "StatQuest: DESeq2, part 1, Library Normalization", https://youtu.be/UFB993xufUU. Problems with library normalization due to library complexity, DESeq2 processing to obtain scaling factors. 

- Video (21 min) "StatQuest: edgeR and DESeq2, part 2 - Independent Filtering", https://youtu.be/Gi0JdrxRq5s. FDR and issues with too many null tests, edgeR filtering by hard cutoff CPM, DESeq2 filtering by quantile CPM plot. Code to calculate DESeq2-like CPM cutoff for edgeR, https://statquest.org/2017/05/16/statquest-filtering-genes-with-low-read-counts/

- Video Naomi Altman, RNA-seq, differential analysis, https://www.youtube.com/watch?v=IWsP-008W0Q&feature=youtu.be

- `Asa_Bjorklund_RNAseqQC_141029.pdf` - RNA-seq QC

- `DiffExp_Estelle_Proux_Vera_201410.pdf` - RNA-seq differential expression, very thorough. 

- `kingsford-regulatory-genomics-salmon.pdf` - Accurate, Fast, and Model-Aware Transcript Expression Quantification with Salmon. https://www.youtube.com/watch?v=TMLIxwDP7sk&list=PLgKuh-lKre10fYXnD-8ghi9b9xl83CKet&index=3


## ToDo 

- Slides from http://combine-australia.github.io/RNAseq-R/

- `RNA-seq_workflow.pdf` - SummarizedExperiment and RNA-seq data manipulation in Bioconductor, https://www.bioconductor.org/help/course-materials/2017/OMRF/B3_RNASeq_Workflow.html


- RNA-seq differential expression with Bioconductor by Davide Risso, https://github.com/drisso/rnaseq_meetup

- DiffSplice, rMATS (shen), EBSeq

- With the TPM normalization, expression values are scaled so that their sum is always 10^6 for each sample. This approach allows transcript proportions to be comparable among samples. However, in case the total mRNA of a sample is dominated by the expression of only a few genes, the remaining fraction of genes will be characterized by especially low expression values. This effect applies only to RNA-seq data, not to microarray data, as RNA-seq does not have an upper limit in its dynamic range (Bullard J.H. Purdom E. Hansen K.D. Dudoit S. Evaluation of statistical methods for normalization and differential expression in mRNA-seq experiments.BMC Bioinformatics. 2010; 11: 94).

### Compositional data analysis

- Fernandes, Andrew D., Jennifer Ns Reid, Jean M. Macklaim, Thomas A. McMurrough, David R. Edgell, and Gregory B. Gloor. “Unifying the Analysis of High-Throughput Sequencing Datasets: Characterizing RNA-Seq, 16S RRNA Gene Sequencing and Selective Growth Experiments by Compositional Data Analysis.” Microbiome 2 (2014): 15. https://doi.org/10.1186/2049-2618-2-15. - ALDEx2 R package - a compositional data analysis tool that uses Bayesian methods to infer technical and statistical errors. Works with RNA-seq, microbiome, and other compositional data. Distinction between absolute counts and compositional data. Counts are converted to probabilities by Monte Carlo sampling (128 by default) from the Dirichlet distribution with a uniform prior. Centered log-ratio transformation, clr - divide by the geometric mean. https://bioconductor.org/packages/release/bioc/html/ALDEx2.html

- Quinn, Thomas P., Ionas Erb, Mark F. Richardson, and Tamsyn M. Crowley. “Understanding Sequencing Data as Compositions: An Outlook and Review.” Bioinformatics (Oxford, England) 34, no. 16 (August 15, 2018): 2870–78. https://doi.org/10.1093/bioinformatics/bty175. - Compositional analysis review as it applies to sequencing data. Introduction with examples, pitfalls of treating compositional data as independent. Normalization to library size as a way to open compositional data for univariate analysis. Aitchinson work, his additive and centered log-ratio transformations as a way to allow tor measurement of Euclidean distances, measurement of association to avoid spurious correlations, principal component analysis. ALDEx2 package implementing those transformations and analysis, using the Dirichlet sampling to avoid zeros. Interpretation difficulties.

- Pitfalls of Data Normalization. https://towardsdatascience.com/pitfalls-of-data-normalization-bf05d65f1f4c

- Van Den Berge, Koen, Katharina Hembach, Charlotte Soneson, Simone Tiberi, Lieven Clement, Michael I Love, Rob Patro, and Mark Robinson. “RNA Sequencing Data: Hitchhiker’s Guide to Expression Analysis.” Preprint. PeerJ Preprints, November 24, 2018. https://doi.org/10.7287/peerj.preprints.27283v2. - RNA-seq review. Experimental technology (sample preparation, Illumina/PacBio/ONT sequencing, experiment design, sample size, sequencing depth). Applications (differential gene expression, transcriptome assembly, allele-specific expression, eQTLs, alternative splicing, pathways, networks). Computation (genome/transcriptome alignment, gene/transcript-level quantification, RSEM). Tools and their brief description. Differential expression strategies (filtering and normalizaition, modeling and estimation, statistical inference, multiple testing, batch effects). Variants of differential analysis (differential transcript/exon usage, tools). Long-read transcript sequencing technologies.

- `TPMcalculator` - converts gene counts to TPM using transcript information from a GTF file. TPM vs. FPKM correlation for validation. C/C++ command line tool, Docker image, CWL workflow. https://github.com/ncbi/TPMCalculator
    - Vera Alvarez, Roberto, Lorinc Sandor Pongor, Leonardo Mariño-Ramírez, and David Landsman. “TPMCalculator: One-Step Software to Quantify MRNA Abundance of Genomic Features.” Edited by Bonnie Berger. Bioinformatics 35, no. 11 (June 1, 2019): 1960–62. https://doi.org/10.1093/bioinformatics/bty896.

- Ribo-depletion in RNA-Seq – Which ribosomal RNA depletion method works best? https://www.rna-seqblog.com/ribo-depletion-in-rna-seq-which-ribosomal-rna-depletion-method-works-best/

- Li, Jun, and Robert Tibshirani. “Finding Consistent Patterns: A Nonparametric Approach for Identifying Differential Expression in RNA-Seq Data.” Statistical Methods in Medical Research 22, no. 5 (October 2013): 519–36. https://doi.org/10.1177/0962280211428386. - Non-parametric statistics for differential expression analysis, including FDR estimation. Compared with Poisson- of negative binomial models, perform better in case of outliers. Lower power. Implemented in samr https://cran.r-project.org/web/packages/samr/index.html and npSeq https://github.com/joey711/npSeq

