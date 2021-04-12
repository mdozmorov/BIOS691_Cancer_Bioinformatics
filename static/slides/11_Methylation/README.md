## `hw` - homework material

## `lab` - in-class lab material

- minfi.Rmd - from https://kasperdanielhansen.github.io/genbioconductor/

- `MethylationAnalysisWorkflow`, https://github.com/Oshlack/MethylationAnalysisWorkflow
- Maksimovic, Jovana, Belinda Phipson, and Alicia Oshlack. “A Cross-Package Bioconductor Workflow for Analysing Methylation Array Data.” F1000Research 5 (2016): 1281. https://doi.org/10.12688/f1000research.8839.3. - Methylation analysis tutorial using multiple packages. Beta and M-values. Differential methylation steps, from data loading, QC, to normalization, filtering, differential CpG and region detection, visualization, functional enrichment analysis, differential variability analysis.https://github.com/Oshlack/MethylationAnalysisWorkflow,  https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html



## References

- Pidsley, R., et. al., and Susan J. Clark. "[Critical Evaluation of the Illumina MethylationEPIC BeadChip Microarray for Whole-Genome DNA Methylation Profiling](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1)." Genome Biology 2016

- Pan D., et. al. "[Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587)." BMC Bioinformatics, 2010

- Bock, Christoph, Eleni M Tomazou, Arie B Brinkman, Fabian Müller, Femke Simmer, Hongcang Gu, Natalie Jäger, Andreas Gnirke, Hendrik G Stunnenberg, and Alexander Meissner. “[Quantitative Comparison of Genome-Wide DNA Methylation Mapping Technologies](https://doi.org/10.1038/nbt.1681).” Nature Biotechnology, (October 2010) - Methylation intro, technology. Software tools, tables. Quality control and problems. Differential analysis. Public repositories.

- Krueger, Felix, Benjamin Kreck, Andre Franke, and Simon R Andrews. “[DNA Methylome Analysis Using Short Bisulfite Sequencing Data](https://www.nature.com/articles/nmeth.1828).” Nature Methods, (January 30, 2012) - Methylation intro, technologies to measure. Alignment problems, QC considerations, processing workflow. Theoretical, references.

- Wreczycka, Katarzyna, Alexander Gosdschan, Dilmurat Yusuf, Bjoern Gruening, Yassen Assenov, and Altuna Akalin. “[Strategies for Analyzing Bisulfite Sequencing Data](https://www.sciencedirect.com/science/article/pii/S0168165617315936),” August 9, 2017 - Bisufite sequencing data analysis steps. Intro into methylation. Refs to packages. 

- Robinson, Mark D., Abdullah Kahraman, Charity W. Law, Helen Lindsay, Malgorzata Nowicka, Lukas M. Weber, and Xiaobei Zhou. “[Statistical Methods for Detecting Differentially Methylated Loci and Regions](https://doi.org/10.3389/fgene.2014.00324).” Frontiers in Genetics (2014) - Methylation review, technology, databases, experimental design, statistics and tools for differential methylation detection, beta-binomial distribution, cell type deconvolution.

- Krueger, Felix, and Simon R. Andrews. “[Bismark: A Flexible Aligner and Methylation Caller for Bisulfite-Seq Applications](https://doi.org/10.1093/bioinformatics/btr167).” Bioinformatics, (June 1, 2011) - Bismark paper. stranded and unstranded BS sequencing. Conversion of reads, genomes, best alignment strategy. https://www.bioinformatics.babraham.ac.uk/projects/bismark/

- Xi, Yuanxin, and Wei Li. “[BSMAP: Whole Genome Bisulfite Sequence MAPping Program](https://doi.org/10.1186/1471-2105-10-232).” BMC Bioinformatics, (2009) - BSMAP paper. Bisulphite conversion technology introduction, problems. BSMAP algorithm. Very good figures explaining all steps. [BSMAP software](https://code.google.com/archive/p/bsmap/)

- Chen, Yunshun, Bhupinder Pal, Jane E. Visvader, and Gordon K. Smyth. “[Differential Methylation Analysis of Reduced Representation Bisulfite Sequencing Experiments Using EdgeR](https://doi.org/10.12688/f1000research.13196.1).” F1000Research (November 28, 2017) - RRBS differential methylation analysis. Methylation intro. R code tutorial. 

- Teschendorff, Andrew E., and Caroline L. Relton. “[Statistical and Integrative System-Level Analysis of DNA Methylation Data](https://doi.org/10.1038/nrg.2017.86).” Nature Reviews Genetics, November 13, 2017 - Deconvolution of methylation profiles. Reference-based, reference-free, semi-reference-free. Table 1 - tools

- Methylation statistics packages: Table 2 in Liu, Hongbo, Song Li, Xinyu Wang, Jiang Zhu, Yanjun Wei, Yihan Wang, Yanhua Wen, et al. “[DNA Methylation Dynamics: Identification and Functional Annotation](https://www.ncbi.nlm.nih.gov/pubmed/27515490).” Briefings in Functional Genomics, 2016

- Mark D. Robinson et al., “[Statistical Methods for Detecting Differentially Methylated Loci and Regions](https://doi.org/10.3389/fgene.2014.00324),” Frontiers in Genetics 5 (2014) - Methylation nethods review. From data, experimental design to software tools for finding differentially methylated regions.

- Shafi, Adib, Cristina Mitrea, Tin Nguyen, and Sorin Draghici. “[A Survey of the Approaches for Identifying Differential Methylation Using Bisulfite Sequencing Data](https://doi.org/10.1093/bib/bbx013).” Briefings in Bioinformatics, March 8, 2017 - Review of differential methylation methods and 22 tools. Categorized by approaches. [Pros and cons of each approach](https://academic.oup.com/view-large/figure/121909042/bbx013f3.tif), [Table 1. Summary of the important characteristics of the 22 surveyed approaches](https://academic.oup.com/view-large/121909044), [Table 2. Comparison of the available implementations of the 22 surveyed approaches](https://academic.oup.com/view-large/121909057)

- Wreczycka, Katarzyna, Alexander Gosdschan, Dilmurat Yusuf, Björn Grüning, Yassen Assenov, and Altuna Akalin. “[Strategies for Analyzing Bisulfite Sequencing Data](https://doi.org/10.1016/j.jbiotec.2017.08.007).” Journal of Biotechnology (November 2017) - Review of bisulfite sequencing technology, data preprocessing, and analysis methods and tools. Differential methylation.

- Venet, D., F. Pecasse, C. Maenhaut, and H. Bersini. “[Separation of Samples into Their Constituents Using Gene Expression Data](https://academic.oup.com/bioinformatics/article/17/suppl_1/S279/262438).” Bioinformatics (2001) - Statistical derivation of deconvolution. 

### Misc

- R annotation ana analysis packages for Illumina methylation arrays, [http://www.hansenlab.org/software.html](http://www.hansenlab.org/software.html)

- Fast and accurante alignment of BS-Seq reads. [https://github.com/brentp/bwa-meth/](https://github.com/brentp/bwa-meth/)

- https://github.com/crazyhottommy/DNA-methylation-analysis - notes on DNA methylation analysis (arrays and sequencing data)


## ToDo - todo list
