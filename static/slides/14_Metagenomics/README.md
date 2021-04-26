- `01_Metagenomics.Rmd` - Overview of metagenomics analysis, from 16S sequencing, OTU clustering to shotgun sequencing and analysis statistics


## `lab` - in-class lab material

- RDP - Ribosomal Database Project. Assignment generator to generate sequences to play with tools. http://rdp.cme.msu.edu/index.jsp

- `/Users/mdozmorov/Documents/Data/Teaching/microbiome_mothur` - data used in mothur paper, Microbial diversity in the deep sea and the underexplored "rare biosphere" (2006); PNAS, 103(32):12115-20. From http://jbpc.mbl.edu/research_supplements/g454/20060412-private/

- `F1000_workflow` - Microbiome workflow. RSV instead of OTU. Data preprocessing from raw reads. DADA2 pipeline, ASV summary tables using RDP (Greengenes and SILVA are available), phylogenetic tree reconstruction (pangorn). phyloseq downstream analysis, from filtering to agglomeration, transformation, various ordination visualizations (from PCoA, DPCoA, rank PCA, to CCA), supervised learning, graph-based visualization and testing, multi-omics analyses. https://github.com/spholmes/F1000_workflow
    - Callahan, Ben J., Kris Sankaran, Julia A. Fukuyama, Paul J. McMurdie, and Susan P. Holmes. “Bioconductor Workflow for Microbiome Data Analysis: From Raw Reads to Community Analyses.” F1000Research 5 (2016): 1492. https://doi.org/10.12688/f1000research.8986.2.


## References

- Morgan, Xochitl C., and Curtis Huttenhower. “[Chapter 12: Human Microbiome Analysis](https://doi.org/10.1371/journal.pcbi.1002808).” PLoS Computational Biology 8, no. 12 (December 27, 2012) - Human Microbiome Analysis review. History, technologies, 16S, OTUs, alpha diversity (Shannon index), beta-diversity (Bray-Curtis dissimilarity) formulas. Includes exercises.

- Goodrich, Julia K., Sara C. Di Rienzi, Angela C. Poole, Omry Koren, William A. Walters, J. Gregory Caporaso, Rob Knight, and Ruth E. Ley. “[Conducting a Microbiome Study](https://doi.org/10.1016/j.cell.2014.06.037).” Cell, (July 17, 2014) - A good review about conducting a microbiome study. Experimental and computational aspects. 16S, ITS, OTU, alpha/beta diversity, Bray-Curtis and UniFrac distance, PCoA.

- Franzosa, Eric A., Tiffany Hsu, Alexandra Sirota-Madi, Afrah Shafquat, Galeb Abu-Ali, Xochitl C. Morgan, and Curtis Huttenhower. “[Sequencing and beyond: Integrating Molecular ‘omics’ for Microbial Community Profiling](https://doi.org/10.1038/nrmicro3451).” Nature Reviews Microbiology, (June 2015) - Metagenomics technology, methods and analysis review. Multi-omics integrative analyses, statistical considerations. Well-written, references.

- Ghodsi, Mohammadreza, Bo Liu, and Mihai Pop. “[DNACLUST: Accurate and Efficient Clustering of Phylogenetic Marker Genes](https://doi.org/10.1186/1471-2105-12-271).” BMC Bioinformatics (June 30, 2011) - DNACLUST - metagenomics clustering of 16S sequencing. Comparison with CD-HIT and UCLUST. Recruit sequences within the radius around cluster seed. Explanation of distance, Needleman-Wunsch algorithm.

- Grice, Elizabeth A., and Julia A. Segre. “[The Human Microbiome: Our Second Genome](https://doi.org/10.1146/annurev-genom-090711-163814).” Annual Review of Genomics and Human Genetics, (September 22, 2012) - Review on microbiome, 16S RNA, tools, human microbiome project.

- Kong, Heidi H. “[Skin Microbiome: Genomics-Based Insights into the Diversity and Role of Skin Microbes](https://doi.org/10.1016/j.molmed.2011.01.013).” Trends in Molecular Medicine, (June 2011) - Skin microbiome review. Introduction about microbiome sequencing (culture and direct), 16S sequencing, other markers. Definitions. Skin-specific findings.

- McDonald, Daniel, Zhenjiang Xu, Embriette R. Hyde, and Rob Knight. “[Ribosomal RNA, the Lens into Life](https://rnajournal.cshlp.org/content/21/4/692).” RNA, (April 2015) - Short review of 16S rRNA history, databases (RDP< Greengenes, SILVA), QIIME

- Huson, D. H., A. F. Auch, J. Qi, and S. C. Schuster. “[MEGAN Analysis of Metagenomic Data](http://www.genome.org/cgi/doi/10.1101/gr.5969107).” Genome Research, (February 6, 2007) - MEGAN paper. Intro into metagenomics, sequencing, alignment analysis. LCA (lowest common ancestor) algorithm.

- Kuczynski, Justin, Christian L. Lauber, William A. Walters, Laura Wegener Parfrey, José C. Clemente, Dirk Gevers, and Rob Knight. “[Experimental and Analytical Tools for Studying the Human Microbiome](https://doi.org/10.1038/nrg3129).” Nature Reviews Genetics, (December 16, 2011)

- Hamady, Micah, and Rob Knight. “[Microbial Community Profiling for Human Microbiome Projects: Tools, Techniques, and Challenges](http://www.genome.org/cgi/doi/10.1101/gr.085464.108).” Genome Research, (July 2009) - Introduction, background of metagenomics. 16S vs. sequencing, experimental questions like read length, sampling depth, how to analyze.

- Schloss, P. D., S. L. Westcott, T. Ryabin, J. R. Hall, M. Hartmann, E. B. Hollister, R. A. Lesniewski, et al. “[Introducing Mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities](https://doi.org/10.1128/AEM.01541-09).” Applied and Environmental Microbiology, (December 1, 2009) - mothur paper. https://www.mothur.org/

- Wang, Q., G. M. Garrity, J. M. Tiedje, and J. R. Cole. “[Naive Bayesian Classifier for Rapid Assignment of RRNA Sequences into the New Bacterial Taxonomy](https://doi.org/10.1128/AEM.00062-07).” Applied and Environmental Microbiology, (August 15, 2007) - Naive Bayes classifier used in Ribosome Database Project (RDP). Good methods statistical description.

- Sczyrba, Alexander, Peter Hofmann, Peter Belmann, David Koslicki, Stefan Janssen, Johannes Dröge, Ivan Gregor, et al. “[Critical Assessment of Metagenome Interpretation-a Benchmark of Metagenomics Software](https://doi.org/10.1038/nmeth.4458).” Nature Methods, (November 2017) - Critical Assessment of Metagenome Interpretation (CAMI) paper. Assessment of microbial genome assemblers, effect of sequencing depth, strain diversity, taxonomic binning, etc., recommendations on software and best practices.

- Letunic, Ivica, and Peer Bork. “[Interactive Tree Of Life (ITOL): An Online Tool for Phylogenetic Tree Display and Annotation](https://doi.org/10.1093/bioinformatics/btl529).” Bioinformatics, (January 1, 2007) - iTOL - manipulation and visualization of phylogenetic trees. New Hampshire and Newick formats, displaying, pruning, annotating.  https://itol.embl.de/

- McMurdie, Paul J., and Susan Holmes. “[Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible](https://doi.org/10.1371/journal.pcbi.1003531).” PLoS Computational Biology, (April 3, 2014) - Rarefying microbiome data or using proportions is wrong, statistical arguments. Variance stabilization (DESeq) and upper-quartile log-fold change normalization (edgeR) perform well to normalize the data. Negative binomial and zero-inflated Gaussian mixture models are recommended to test for differential abundance, differential abundance test in metagenomeSeq package also performs well. Rarefying lead to high proportion of false positives. Importance of filtering.

- [Microbiome Discovery](https://www.youtube.com/playlist?list=PLOPiWVjg6aTzsA53N19YqJQeZpSCH9QPc) - video course playlist about microbiome, from biology to statistics, by Dan Knights.

### Metagenomic resources

- [awesome-microbes](https://github.com/stevetsa/awesome-microbes) - List of software packages (and the people developing these methods) for microbiome (16S), metagenomics (WGS, Shot-gun sequencing), and pathogen identification/detection/characterization. 

- Schiffer, Lucas, Rimsha Azhar, Lori Shepherd, Marcel Ramos, Ludwig Geistlinger, Curtis Huttenhower, Jennifer B Dowd, Nicola Segata, and Levi Waldron. “[HMP16SData: Efficient Access to the Human Microbiome Project through Bioconductor](https://doi.org/10.1101/299115),” August 29, 2018 - [HMP16SData](https://bioconductor.org/packages/HMP16SData/) - SummarizedExperiment of 16S sequencing data (counts) for V13 and V35 variable regions with clinical annotations (visit number, sex, run center, body site, and body subsite) and an option to attach controlled access clinical annotations. Compatible with phylosec. 

- [Microbiome Helper](https://github.com/LangilleLab/microbiome_helper/wiki) - wrapper scripts and tutorials for metagenomics analysis. Main paper: Comeau, André M., Gavin M. Douglas, and Morgan G. I. Langille. “[Microbiome Helper: A Custom and Streamlined Workflow for Microbiome Research](https://doi.org/10.1128/mSystems.00127-16).” MSystems, (February 28, 2017)


## ToDo

- Bray-Curtis dissimilarity metric, https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity

- Microbiome data analysis by the Waldron lab, https://github.com/waldronlab/MicrobiomeWorkshop/tree/BiocNYC-2017-12-15


## `misc`

- `Microbiome module` - Source: https://github.com/TheJacksonLaboratory/JAXBD2K-ShortCourse

- `10-A-Kumar_slide_show.pdf` - https://www.uab.edu/medicine/camac/images/10-A-Kumar_slide_show.pdf
