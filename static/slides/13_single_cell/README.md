## `hw` - homework material

## `lab` - in-class lab material

https://github.com/BIMSBbioinfo/SingleCell_2018

## References

### Technology

- Svensson, Valentine, Roser Vento-Tormo, and Sarah A Teichmann. “[Exponential Scaling of Single-Cell RNA-Seq in the Past Decade](https://doi.org/10.1038/nprot.2017.149).” Nature Protocols, (March 1, 2018) - scRNA-seq technology development.

- Wang, Yong, and Nicholas E. Navin. “[Advances and Applications of Single-Cell Sequencing Technologies](https://doi.org/10.1016/j.molcel.2015.05.005).” Molecular Cell, (May 21, 2015) - Single cell sequencing review, all technologies. Whole genome amplification.

- Zheng, Grace X. Y., Jessica M. Terry, Phillip Belgrader, Paul Ryvkin, Zachary W. Bent, Ryan Wilson, Solongo B. Ziraldo, et al. “[Massively Parallel Digital Transcriptional Profiling of Single Cells](https://doi.org/10.1038/ncomms14049).” Nature Communications 8 (January 16, 2017) - 10X technology. Details of each wet-lab step, sequencing, and basic computational analysis.

- Dal Molin, Alessandra, and Barbara Di Camillo. “[How to Design a Single-Cell RNA-Sequencing Experiment: Pitfalls, Challenges and Perspectives](https://doi.org/10.1093/bib/bby007).” Briefings in Bioinformatics, January 31, 2018 - single-cell RNA-seq review. From experimental to bioinformatics steps. Cell isolation (FACS, microfluidics, droplet-based), mRNA capture, RT and amplification (poly-A tails, template switching, IVT), quantitative standards (spike-ins, UMIs), transcript quantification, normalization, batch effect removal, visualization.

### Methods

- Camara, Pablo G. “[Methods and Challenges in the Analysis of Single-Cell RNA-Sequencing Data](https://doi.org/10.1016/j.coisb.2017.12.007).” Current Opinion in Systems Biology (February 2018) - scRNA-seq conscise review of computational analysis steps.

- Bacher, Rhonda, and Christina Kendziorski. “[Design and Computational Analysis of Single-Cell RNA-Sequencing Experiments](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y).” Genome Biology (April 7, 2016) - scRNA-seq analysis. Table 1 - categorized tools and description. Normalization, noise reduction, clustering, differential expression, pseudotime ordering. 

- Vallejos, Catalina A, Davide Risso, Antonio Scialdone, Sandrine Dudoit, and John C Marioni. “[Normalizing Single-Cell RNA Sequencing Data: Challenges and Opportunities](https://doi.org/10.1038/nmeth.4292).” Nature Methods, (May 15, 2017) - single-cell RNA-seq normalization methods. Noise, sparsity. Cell- and gene-specific effects. Bulk RNA-seq normalization methods don't work well. Overview of RPKM, TPM, TMM, DESeq normalizations. Spike-in based normalization methods.

- Freytag, Saskia, Luyi Tian, Ingrid Lönnstedt, Milica Ng, and Melanie Bahlo. “[Comparison of Clustering Tools in R for Medium-Sized 10x Genomics Single-Cell RNA-Sequencing Data](https://doi.org/10.12688/f1000research.15809.1).” F1000Research (August 15, 2018) - Comparison of 12 scRNA-seq clustering tools on medium-size 10X genomics datasets. Adjusted Rand index, normalized mutual information, homogeneity score, stability by random sampling of samples or genes. Seurat and Cell Ranger outperform other methods.

### Statistics

- Brennecke, Philip, Simon Anders, Jong Kyoung Kim, Aleksandra A Kołodziejczyk, Xiuwei Zhang, Valentina Proserpio, Bianka Baying, et al. “[Accounting for Technical Noise in Single-Cell RNA-Seq Experiments](https://doi.org/10.1038/nmeth.2645).” Nature Methods, (September 22, 2013) - Single-cell noise. Technical, biological. Use spike-ins to estimate noise. Can be approximated with Poisson distribution.

- Grün, Dominic, Lennart Kester, and Alexander van Oudenaarden. “[Validation of Noise Models for Single-Cell Transcriptomics](https://doi.org/10.1038/nmeth.2930).” Nature Methods, (June 2014) - Quantification of sampling noise and global cell-to-cell variation in sequencing efficiency. Three noise models. 4-bases-long UMIs are sufficient for transcript quantification, improve CV. Negative binomial component of expressed genes. Statistics of transcript counting from UMIs, negative binomial distribution, noise models

- Risso, Davide, Fanny Perraudeau, Svetlana Gribkova, Sandrine Dudoit, and Jean-Philippe Vert. “[ZINB-WaVE: A General and Flexible Method for Signal Extraction from Single-Cell RNA-Seq Data](https://doi.org/10.1101/125112).” BioRxiv, January 1, 2017 - Zero-inflated negative binomial model for normalization, batch removal, and dimensionality reduction. Extends the RUV model with more careful definition of "unwanted" variation as it may be biological. Good statistical derivations in Methods. Refs to real and simulated scRNA-seq datasets

- Ding, Bo, Lina Zheng, Yun Zhu, Nan Li, Haiyang Jia, Rizi Ai, Andre Wildberg, and Wei Wang. “[Normalization and Noise Reduction for Single Cell RNA-Seq Experiments](https://doi.org/10.1093/bioinformatics/btv122).” Bioinformatics, (July 1, 2015) - Fitting gamma distribution to log2 read counts of known spike-in ERCC controls, predicting RNA concentration from it.

- Bacher, Rhonda, Li-Fang Chu, Ning Leng, Audrey P Gasch, James A Thomson, Ron M Stewart, Michael Newton, and Christina Kendziorski. “[SCnorm: Robust Normalization of Single-Cell RNA-Seq Data](https://doi.org/10.1038/nmeth.4263).” Nature Methods, (April 17, 2017) - SCnorm - normalization for single-cell data. Quantile regression to estimate the dependence of transcript expression on sequencing depth for every gene. Genes with similar dependence are then grouped, and a second quantile regression is used to estimate scale factors within each group. Within-group adjustment for sequencing depth is then performed using the estimated scale factors to provide normalized estimates of expression. Good statistical methods description

- Finak, Greg, Andrew McDavid, Masanao Yajima, Jingyuan Deng, Vivian Gersuk, Alex K. Shalek, Chloe K. Slichter, et al. “[MAST: A Flexible Statistical Framework for Assessing Transcriptional Changes and Characterizing Heterogeneity in Single-Cell RNA Sequencing Data](https://doi.org/10.1186/s13059-015-0844-5).” Genome Biology (December 10, 2015) - MAST, scRNA-seq DEG analysis. CDR - the fraction of genes that are detectably expressed in each cell - added to the hurdle model that explicitly parameterizes distributions of expressed and non-expressed genes.


### Workflows and tools

- Amezquita, Robert A., Vincent J. Carey, Lindsay N. Carpp, Ludwig Geistlinger, Aaron TL Lun, Federico Marini, Kevin Rue-Albrecht, et al. “[Orchestrating Single-Cell Analysis with Bioconductor](https://doi.org/10.1038/s41592-019-0654-x).” Nature Methods, 13 September 2019 - scRNA-seq analysis overview within Bioconductor ecosystem. SingleCellexperiment, scran and scater examples. Table S1 - summary of packages for data input, infrastructure, QC, integration, dimensionality reduction, clustering, pseudotime, differential expression, functional enrichment, simulation.

- Poirion, Olivier B., Xun Zhu, Travers Ching, and Lana Garmire. “[Single-Cell Transcriptomics Bioinformatics and Computational Challenges](https://doi.org/10.3389/fgene.2016.00163).” Frontiers in Genetics (2016) - single-cell RNA-seq review. Workflow, table 1 - tools, table 2 - data. Workflow sections describe each tools. References to all tools.

- Lun, Aaron T. L., Davis J. McCarthy, and John C. Marioni. “[A Step-by-Step Workflow for Low-Level Analysis of Single-Cell RNA-Seq Data with Bioconductor](https://doi.org/10.12688/f1000research.9501.2).” F1000Research (2016) - scRNA-seq analysis, from count matrix. Noise. Use of ERCCs and UMIs. QC metrics - library size, expression level, mitochondrial genes. Tests for batch effect. Finding highly variable genes, clustering. Several examples, with R code. [A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor](https://bioconductor.org/help/workflows/simpleSingleCell/), [Methods for Single-Cell RNA-Seq Data Analysis](https://bioconductor.org/packages/devel/bioc/html/scran.html)

- Perraudeau, Fanny, Davide Risso, Kelly Street, Elizabeth Purdom, and Sandrine Dudoit. “[Bioconductor Workflow for Single-Cell RNA Sequencing: Normalization, Dimensionality Reduction, Clustering, and Lineage Inference](https://doi.org/10.12688/f1000research.12122).” F1000Research (July 21, 2017) - single-cell RNA-seq pipeline. Post-processing analysis in R, uzinng SummarizedExperiment object, scone for QC, ZINB-WaVE for normalization and dim. reduction, clusterExperiment for clustering, slingshot for temporal inference, differential expression using generalized additive model

- McCarthy, Davis J., Kieran R. Campbell, Aaron T. L. Lun, and Quin F. Wills. “[Scater: Pre-Processing, Quality Control, Normalization and Visualization of Single-Cell RNA-Seq Data in R](https://doi.org/10.1093/bioinformatics/btw777).” Bioinformatics, (April 15, 2017)

- Zappia, Luke, Belinda Phipson, and Alicia Oshlack. “[Exploring the Single-Cell RNA-Seq Analysis Landscape with the ScRNA-Tools Database](https://doi.org/10.1101/206573).” BioRxiv, January 1, 2018 - [www.scrna-tools.org](www.scrna-tools.org) - structured collection of scRNA-seq tools, from visualization, specialized tools, to full pipelines

## scRNA-seq resources

- [scRNA-seq analysis notes by Ming Tang, tutorials, tools, papers, code snippets](https://github.com/crazyhottommy/scRNAseq-analysis-notes)

- [List of software packages for multi-omics analysis by Mike Love](https://github.com/mikelove/awesome-multi-omics)

- [List of software packages for single-cell data analysis, including RNA-seq, ATAC-seq, etc by Sean Davis](https://github.com/seandavi/awesome-single-cell)

- [Analysis of single cell RNA-seq data course, Cambridge University, UK](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html). [Videos](https://www.youtube.com/playlist?list=PLEyKDyF1qdOYAhwU71qlrOXYsYHtyIu8n)

- [Single cell data portal](https://portals.broadinstitute.org/single_cell)

- [Conquer DB of scRNA-seq datasets](http://imlspenticton.uzh.ch:3838/conquer/)



## ToDo - todo list

- A comprehensive list of single-cell multi-omics methods, software, companies. https://docs.google.com/spreadsheets/d/1IPe2ozb1Mny8sLvJaSE57RJr3oruiBoSudAVhSH-O8M/edit#gid=0

- Figure depicting the breadth of multimodal scRNA-seq technologies. https://github.com/arnavm/multimodal-scRNA-seq

- `scSeq.pdf` - Single-cell sequencing, http://www.haowulab.org/teaching/bioc/bioc.html

- `single-cell-pseudotime` - Single-cell RNA-seq pseudotime estimation algorithms, comprehensive collection of links to software and accompanying papers. https://github.com/agitter/single-cell-pseudotime - to read more

- "Getting started with Monocle", https://davetang.org/muse/2017/10/01/getting-started-monocle/

- powsimR: Power analysis for bulk and single cell RNA-seq experiments, https://github.com/bvieth/powsimR

- `lab/SingleCell_2018` - Repository for Single Cell Analysis Course in MDC

## References

- Single-cell ATAC-seq: Buenrostro, Jason D., Beijing Wu, Ulrike M. Litzenburger, Dave Ruff, Michael L. Gonzales, Michael P. Snyder, Howard Y. Chang, and William J. Greenleaf. “Single-Cell Chromatin Accessibility Reveals Principles of Regulatory Variation.” Nature 523, no. 7561 (July 23, 2015): 486–90. https://doi.org/10.1038/nature14590.

- Hwang, Byungjin, Ji Hyun Lee, and Duhee Bang. “Single-Cell RNA Sequencing Technologies and Bioinformatics Pipelines.” Experimental & Molecular Medicine 50, no. 8 (August 2018). https://doi.org/10.1038/s12276-018-0071-8. - scRNA-seq technology and computational analysis review. Single cell isolation techniques, library preparation methods, computational analyses, R/FPKM/TPM, batch and cell cycle effect, cell type identification, clustering, network inference, trajectory inference. Figures for teaching.

- scRNA-seq technologies, including multi-omics, from the Discussion at https://www.nature.com/articles/s41467-018-07771-0
    - In recent years, other methods, such as DNase-seq32, MNase-seq33 and NOMe-seq34,35, have investigated chromatin status at the single cell level. However, due to its simplicity and reliability, ATAC-seq currently remains the most popular technique for chromatin profiling. Several recent studies have demonstrated the power of using scATAC-seq for investigating regulatory principles, e.g. brain development4,9, Mouse sci-ATAC-seq Atlas36 and pseudotime inference37. The combined multi-omics approaches also began to emerge, such as sci-CAR-seq38, scCAT-seq39 and piATAC-seq8.

- A workflow for single cell RNA-seq data analysis, http://research.fhcrc.org/content/dam/stripe/sun/software/scRNAseq/scRNAseq.html

- Introduction to single-cell RNA-seq technologies, presentation by Lior Pachter. Key figures, references, statistics. Slides, https://figshare.com/articles/Introduction_to_single-cell_RNA-seq_technologies/7704659/1, and notes, https://liorpachter.wordpress.com/2019/02/19/introduction-to-single-cell-rna-seq-technologies/

- Luecken, Malte D., and Fabian J. Theis. “Current Best Practices in Single-Cell RNA-Seq Analysis: A Tutorial.” Molecular Systems Biology 15, no. 6 (June 19, 2019): e8746. https://doi.org/10.15252/msb.20188746. - All steps in scRNA-seq analysis, QC (count depth, number of genes, % mitochondrial), normalization (global, downsampling, nonlinear), data correction (batch, denoising, imputation), feature selection, dimensionality reduction (PCA, diffusion maps, tSNE, UMAP), visualization, clustering (k-means, graph/community detection), annotation, trajectory inference (PAGA, Monocle), differential analysis (DESeq2, EdgeR, MAST), gene regulatory networks. Description of the bigger picture at each step, latest tools, their brief description, references. R-based Scater as the full pipeline for QC and preprocessing, Seurat for downstream analysis, scanpy Python pipeline. Links and refs for tutorials. https://github.com/theislab/single-cell-tutorial

# ToDo

## Hwang, Byungjin, Ji Hyun Lee, and Duhee Bang. “Single-Cell RNA Sequencing Technologies and Bioinformatics Pipelines.” Experimental & Molecular Medicine 50, no. 8 (07 2018): 96. https://doi.org/10.1038/s12276-018-0071-8. - scRNA-seq technology and analysis overview. Isolation techniques (Micromanipulation, FACS, microfluidics, microdroplet). Only 10-20% of transcripts get reverse-transcribed. UMIs, RPKM/FPKM, TPM, TMM, DESeq normalization. Dimensionality reduction. Network generation. Cell hierarchy/pseudotime reconstruction.

The sequencing an entire transcriptome at the level of a single-cell was pioneered by James Eberwine et al.6 and Iscove and colleagues7, who expanded the complementary DNAs (cDNAs) of an individual cell using linear ampli- fication by in vitro transcription and exponential ampli- fication by PCR, respectively.
6. Eberwine, J. et al. Analysis of gene expression in single live neurons. Proc.
Natl. Acad. Sci. USA 89, 3010–3014 (1992).
7. Brady, G., Barbara, M. & Iscove, N. N. Representative in vitro cDNA amplification from individual hemopoietic cells and colonies. Methods Mol. Cell Biol.
2, 17–25 (1990).

Isolation techniques (Micromanipulation, FACS, microfluidics, microdroplet)

It has been well established that due to Poisson sampling, only 10–20% of transcripts will be reverse transcribed at this stage32.
 Islam, S. et al. Quantitative single-cell RNA-seq with unique molecular
identifiers. Nat. Methods 11, 163–166 (2014)


quantification.png Fig. 3 Methods for the quantification of expression in scRNA-seq. a Reads per kilobase (RPK) is defined by multiplying the read counts of an isoform (i) by 1000 and dividing by isoform length. Reads per kilobase per million (RPKM) is defined to compare experiments or different samples (cells) so that additional normalization by the total fragment count is integrated in the denominator term, which is expressed in millions. b The metric TPM takes other isoforms into account, which contrasts with the RPKM metric. This metric quantifies the abundance of isoforms (i) using the RPK fraction across isoforms. c A schematic example illustrates the difference between the RPKM and TPM measures. TPM is efficient for measuring relative abundance because total normalized reads are constant across different cells. d However, we should be careful to interpret the fact that differentially expressed genes can be falsely annotated as a result of overexpression of the other isoforms
 
To overcome the inherent problems in within-sample normalization methods, alternative approaches have been developed52–54. The trimmed mean of M-values (TMM) method and DESeq are the two most popular choices for between-sample normalization. The basic idea behind these frameworks is that highly variable genes dominate the counts, thus skewing the relative abundance in expression profiles. First, TMM picks reference samples, and the other samples are considered test samples. M- values for each gene are calculated as the genes’ log expression ratios between tests to the reference sample. Then, after excluding the genes with extreme M-values, the weighted average of these M-values is set for each test sample. Similar to TMM, DESeq calculates the scaling factor as the median of the ratios of each gene’s read count in the particular sample over its geometric mean across all samples. However, both approaches (TMM, DESeq) will perform poorly when a large number of zero counts are present. A normalization method based on pooling expression values55 were developed to avoid sto- chastic zero counts which is robust to differentially expressed genes in the data. 
52. Robinson, M. D. & Oshlack, A. A scaling normalization method for differential
expression analysis of RNA-seq data. Genome Biol. 11, R25–R25 (2010).
53. Anders, S. & Huber, W. Differential expression analysis for sequence count
data. Genome Biol. 11, R106–R106 (2010).
54. Li, J., Witten, D. M., Johnstone, I. M. & Tibshirani, R. Normalization, testing, and
false discovery rate estimation for RNA-sequencing data. Biostatistics 13,
523–538 (2012).
55. Lun, A. T., Bach, K. & Marioni, J. C. Pooling across cells to normalize single-cell
RNA sequencing data with many zero counts. Genome Biol. 17, 75 (2016).


The concept of “pseudotime” was introduced in the Monocle16 algorithm, which measures a cell’s biological progression (Fig. 5). Here, the notion of “pseudotime” is different from “real time” because cells are sampled all at once. Maximum parsimony is the basic principle that infers cellular dynamics and has been widely used in phylogenetic tree reconstruction in evolutionary biol- ogy86,87. Monocle initially builds graphs in which the nodes represent cells and the edges correspond to each pair of cells. The edge weights are calculated based on the distance between cells in the matrix obtained from dimensionality reduction using independent component analysis (ICA). The minimum spanning tree (MST) algorithm is then applied to search for the longest back- bone. The main limitation of these methods is that the constructed tree is highly complex, and therefore, the user must specify k branches to search.

- `misc/Single-cell_bioinformatics.pdf` - 'Overview of single cell bioinformatics' by Peter Hickley, https://tinyurl.com/scRNA-bioinf

- Multi-sample multi-group scRNA-seq analysis tools - https://github.com/HelenaLC/muscat

- Harvard FAS informatics nanocourse: scRNAseq workshop, by Tommy Tang. https://github.com/crazyhottommy/scRNA-seq-workshop-Fall-2019

- Seurat_object_explained.pptx - https://osf.io/49q2u/, Crazyhottommy

- Reduce Dimensions for Single Cell, https://towardsdatascience.com/reduce-dimensions-for-single-cell-4224778a2d67

- Comparing normalization strategies for single cell genomics, https://medium.com/@nikolay.oskolkov/how-to-normalize-single-cell-a438281ea654

- Tritschler, Sophie, Maren Büttner, David S. Fischer, Marius Lange, Volker Bergen, Heiko Lickert, and Fabian J. Theis. “Concepts and Limitations for Learning Developmental Trajectories from Single Cell Genomics.” Development 146, no. 12 (June 15, 2019): dev170506. https://doi.org/10.1242/dev.170506. - Review of pseudotime/trajectory inference methods for scRNA-seq. Basic concepts (Boxes). Reference to https://dynverse.org/ for methods comparison. Lineage trees, graph abstraction. Limitations, confounding factors and methods to account for them.

- How to Batch Correct Single Cell, https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1

- Single-cell RNA seq module, slides and exercises. https://kkorthauer.org/fungeno2019/

- How Exactly UMAP Works, And why exactly it is better than tSNE, by Nikolay Oskolkov, https://towardsdatascience.com/how-exactly-umap-works-13e3040e1668

- Understanding UMAP, https://pair-code.github.io/understanding-umap/

- Wang, Tianyu, Boyang Li, Craig E. Nelson, and Sheida Nabavi. “Comparative Analysis of Differential Gene Expression Analysis Tools for Single-Cell RNA Sequencing Data.” BMC Bioinformatics 20, no. 1 (December 2019): 40. https://doi.org/10.1186/s12859-019-2599-6. - Comparison of 11 differential gene expression detection methods for scRNA-seq data. Variable performance, poor overlap. Brief description of the statistics of each method. Bulk RNA-sec methods perform well, edgeR is good and fast. [Table 1. Software tools for identifying DE genes using scRNA-seq data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2599-6/tables/1)

- Distributions to model scRNA-seq data:
    - negative binomial distribution, https://divingintogeneticsandgenomics.rbind.io/post/negative-binomial-distribution-in-scrnaseq/
    - negative bionomial distribution in (single-cell) RNAseq, https://divingintogeneticsandgenomics.rbind.io/post/negative-bionomial-distribution-in-single-cell-rnaseq/
    - Modeling single cell RNAseq data with multinomial distribution, https://divingintogeneticsandgenomics.rbind.io/post/modeling-single-cell-rnaseq-data-with-multinomial-distribution/

- Recommendations to properly use t-SNE on large omics datasets (scRNA-seq in particular) to preserve global geometry. Overview of t-SNE, PCA, MDS, UMAP, their similarities, differences, strengths and weaknesses. PCA initialization (first two components are OK), a high learning rate of n/12, and multi-scale similarity kernels. For very large data, increase exagerration. Strategies to align new points on an existing t-SNE plot, aligning two t-SNE visualizations. Extremely fast implementation is FIt-SNE, https://github.com/KlugerLab/FIt-SNE. Code to illustrate the use of t-SNE: https://github.com/berenslab/rna-seq-tsne
    - Kobak, Dmitry, and Philipp Berens. “The Art of Using T-SNE for Single-Cell Transcriptomics.” Nature Communications 10, no. 1 (December 2019): 5416. https://doi.org/10.1038/s41467-019-13056-x.

- Downstream analysis of scRNAseq data. Gene Set enrichment (GSEA) analysis, https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html

- https://ptaus.netlify.com/post/01_tsne_explanation/tsne_explanation/

- Logistic regression for differential expression in scRNA-seq. Works on pseudocounts, TPM, RSEM. Outperforms Monocle, DESeq2, MAST, SCDE. Implementation - R glm function. https://github.com/pachterlab/NYMP_2018
    - Ntranos, Vasilis, Lynn Yi, Páll Melsted, and Lior Pachter. “Identification of Transcriptional Signatures for Cell Types from Single-Cell RNA-Seq.” Preprint. Bioinformatics, February 1, 2018. https://doi.org/10.1101/258566.
