- `07a_Intro.Rmd` - Introduction into differential expression analysis
- `lab/DEG_ttest.R` - Basic t-test and volcano plot
- `lab/volcano.R` - Demo of pretty volcano plot
- `Power_analysis.Rmd` - Power analysis handouts
- `07b_design_matrices.Rmd` - Refresher on design matrices
- `07b_limma.Rmd` - Linear models for microarray data analysis, limma, empirical Bayes 
- `lab/DiffExpr_Limma.Rmd` - limma demo. Uses data from `/Users/mdozmorov/Documents/Work/Teaching/ci-workshop/data/Su_CELs`
- `07c_SAM.Rmd` - SAM, significance analysis of microarrays
- `lab/DiffExpr_SAM.Rmd` - SAM demo, uses GDS858 data. BiomaRt demo at the end.
- `07d_multipletesting.Rmd` - Multiple testing correction
- `lab/Filtering.Rmd` - filtering gene set
- `lab/multitest.R` - multiple testing correction and filtering. Uses `Methyl.RData`
- `07e_batch.Rmd` - Batch effect correction, ComBat, SVA
- `lab/batch_sva_combat.R` - "Batch effects and confounders", by Jeff Leek, [Source](http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html)
- `lab/adjusting_with_linear_models.Rmd` - another batch removal demo, [Source](https://genomicsclass.github.io/book/pages/adjusting_with_linear_models.html), with added `limma::removeBatchEffect`, see [paper](https://www.ncbi.nlm.nih.gov/pubmed/17206142) for data details

## `hw` - homework material

## `lab` - in-class lab material

- `heatmaps.R` - demo of scatterplots and heatmaps
- `Multiple_Testing.R` - Mark Reimers lab
- `permutations.R` - permutation test


## `read` - reading assignment

- Krzywinski, Martin, and Naomi Altman. “[Points of Significance: Comparing Samples—part I](https://doi.org/10.1038/nmeth.2858).” Nature Methods, (March 2014)

- Krzywinski, Martin, and Naomi Altman. “[Points of Significance: Analysis of Variance and Blocking](https://doi.org/10.1038/nmeth.3005).” Nature Methods, (July 2014)

- Krzywinski, Martin, and Naomi Altman. “[Points of Significance: Power and Sample Size](https://doi.org/10.1038/nmeth.2738).” Nature Methods, (November 26, 2013)

- Tong, Tiejun, and Hongyu Zhao. “[Practical Guidelines for Assessing Power and False Discovery Rate for a Fixed Sample Size in Microarray Experiments](https://doi.org/10.1002/sim.3237).” Statistics in Medicine, (May 20, 2008) - Power analysis. t-statistics, FDR types and definitions, then derivation of power calculations. 

- Chapter 4 "Matrix Algebra", Chapter 5 "Linear Models", from [PH525x series - Biomedical Data Science](https://genomicsclass.github.io/book/) course

- [ANOVA, F-test explanation](https://youtu.be/DFJn0O0-AfU) 10m video
- [Design matrices](https://youtu.be/2UYx-qjJGSs) 14m video

- Jelle J. Goeman and Aldo Solari, “[Multiple Hypothesis Testing in Genomics](https://doi.org/10.1002/sim.6082),” Statistics in Medicine, (May 20, 2014)

- [Statistics for Genomics: Advanced Differential Expression](https://youtu.be/QINX3cI7qgk) 24m video by Rafael Irizarry - limma statistics, Bayesian intro

- Goeman, Jelle J., and Aldo Solari. “[Multiple Hypothesis Testing in Genomics](https://doi.org/10.1002/sim.6082).” Statistics in Medicine, (May 20, 2014) - multiple testing review

- Benjamini, Yoav, and Yosef Hochberg. “[Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing](https://www.jstor.org/stable/2346101).” Journal of the Royal Statistical Society, - FDR paper

- Storey, John D., and Robert Tibshirani. “[Statistical Significance for Genomewide Studies](https://doi.org/10.1073/pnas.1530509100).” Proceedings of the National Academy of Sciences of the United States of America, (August 5, 2003) - q-value paper

- Krzywinski, Martin, and Naomi Altman. “[Points of Significance: Comparing Samples—part II](https://doi.org/10.1038/nmeth.2900).” Nature Methods, (March 28, 2014)

- Goh, Wilson Wen Bin, and Limsoon Wong. “Dealing with Confounders in Omics Analysis.” Trends in Biotechnology 36, no. 5 (May 2018): 488–98. https://doi.org/10.1016/j.tibtech.2018.01.013. - Batch effect, p-value issues, conforunding effects, bias, PCA, problems with chi-square enrichment analysis, permutation test pluses and minuses.

- Lazar, Cosmin, Stijn Meganck, Jonatan Taminau, David Steenhoff, Alain Coletta, Colin Molter, David Y. Weiss-Solís, Robin Duque, Hugues Bersini, and Ann Nowé. “Batch Effect Removal Methods for Microarray Gene Expression Data Integration: A Survey.” Briefings in Bioinformatics 14, no. 4 (July 2013): 469–90. https://doi.org/10.1093/bib/bbs037. - Review of batch correction methods, definitions, methodologies

- Leek, Jeffrey T., Robert B. Scharpf, Héctor Corrada Bravo, David Simcha, Benjamin Langmead, W. Evan Johnson, Donald Geman, Keith Baggerly, and Rafael A. Irizarry. “Tackling the Widespread and Critical Impact of Batch Effects in High-Throughput Data.” Nature Reviews. Genetics 11, no. 10 (October 2010): 733–39. https://doi.org/10.1038/nrg2825. - Batch effect, types, sources, examples of wrong conclusions, SVA and ComBat methods

- Johnson, W. Evan, Cheng Li, and Ariel Rabinovic. “Adjusting Batch Effects in Microarray Expression Data Using Empirical Bayes Methods.” Biostatistics (Oxford, England) 8, no. 1 (January 2007): 118–27. https://doi.org/10.1093/biostatistics/kxj037. - ComBat paper. Batch effect removal using Empirical Bayes method.

- Leek, Jeffrey T., and John D. Storey. “Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis.” PLoS Genetics 3, no. 9 (September 2007): 1724–35. https://doi.org/10.1371/journal.pgen.0030161. - SVA paper

- Chen, Chao, Kay Grennan, Judith Badner, Dandan Zhang, Elliot Gershon, Li Jin, and Chunyu Liu. “Removing Batch Effects in Analysis of Expression Microarray Data: An Evaluation of Six Batch Adjustment Methods.” PloS One 6, no. 2 (2011): e17238. https://doi.org/10.1371/journal.pone.0017238. - Comparing batch effect removal software. ComBat is best. Importance of standartization




## `misc` - misc presentations and materials

- `lec10-linear-models-basic-intro` - Full explanation of linear models. Megative binomial https://github.com/Bioconductor/CSAMA/tree/2017/lecture/3-wednesday/lec10-linear-models-basic-intro  
    - `Waldron_linearmodels` - lecture

- `4b_LinearModelsAndEmpiricalBayes.pdf` - "Linear Models and Empirical Bayes Methods for Microarray Data Analysis", Alex Sánchez. http://www.ub.edu/stat/docencia/bioinformatica/microarrays/ADM/slides/4b_LinearModelsAndEmpiricalBayes.pdf

- `Lecture5-Feb11-08.pdf` - SAM, important limma notes by Gordon Smyth. http://odin.mdacc.tmc.edu/~kim/TeachBioinf/AdvStatGE-Prot.htm


- `batch.pdf` - Batch effect, SVA, SVD, PCA. Kasper Hansen, http://www.biostat.jhsph.edu/~khansen/teaching/2014/140.668/batch.pdf

- `cbmb-shortcourse-lecture1.ppt` - full lecture on microarray preprocessing, differential analysis, limma, sam. http://www.biostat.ucsf.edu/cbmb/teaching.html

- `de140688_2014.pdf` - full lecture on differential analysis, MA-plot, GSEA. differential expression, Rafa Irizarry. http://www.biostat.jhsph.edu/~khansen/teaching/2014/140.668/de140688_2014.pdf

- `DifferentialExpressionII.pdf` - Bayesian approach, limma, sam, multiple testing

- `jsm_2016.pdf` - latent variable analysis



- `Satterthwaite.pdf` - derivations. http://apcentral.collegeboard.com/apc/public/repository/ap05_stats_allwood_fin4prod.pdf

- `t.perm.R` - (P. Legendre): t-test for independent samples with permutation test. [Source](http://adn.biol.umontreal.ca/~numericalecology/Rcode/)
- `t.paired.perm.R` - (P. Legendre & G. Blanchet): t-test for related samples with permutation test.
- `Sidak.R` - (P. Legendre): Bonferroni and Sidak corrections for multiple testing.

- `07b_design_matrices.Rmd` - Refresher on design matrices. Done from https://www.youtube.com/watch?v=nk2CQITm_eo&list=PLblh5JKOoLUIcdlgu78MnlATeyx4cEVeR&index=15, https://www.youtube.com/watch?v=NF5_btOaCig&list=PLblh5JKOoLUIcdlgu78MnlATeyx4cEVeR&index=26,  https://www.youtube.com/watch?v=2UYx-qjJGSs&list=PLblh5JKOoLUIcdlgu78MnlATeyx4cEVeR&index=33

## ToDo

- Festing, Michael FW. “On Determining Sample Size in Experiments Involving Laboratory Animals.” Laboratory Animals, n.d., 10. - Sample size, power analysis for two-group comparison in animal experiments. Basics of power analysis, each component explained. Table 1 lists sample size for a given standardized effect size.

Most existing methods for batch correction are based on linear regression. The limma package provides the removeBatchEffect function, which fits a linear model containing a blocking term for the batch structure to the expression values for each gene. Subsequently, the coefficient for each blocking term is set to zero, and the expression values are computed from the remaining terms and residuals, thus yielding a new expression matrix without batch effects. The ComBat method uses a similar strategy but performs an additional step involving empirical Bayes shrinkage of the blocking coefficient estimates. This procedure stabilizes the estimates in the presence of limited replicates by sharing information across genes. Other methods, such as RUVseq and svaseq, are also frequently used for batch correction, but their focus is primarily on identifying unknown factors of variation, for example, those due to unrecorded experimental differences in cell processing. After these factors are identified, their effects can be regressed out as described previously. [from Haghverdi, Laleh, Aaron T L Lun, Michael D Morgan, and John C Marioni. “Batch Effects in Single-Cell RNA-Sequencing Data Are Corrected by Matching Mutual Nearest Neighbors.” Nature Biotechnology, April 2, 2018. https://doi.org/10.1038/nbt.4091]

- Understanding p value, multiple comparisons, FDR and q value, https://divingintogeneticsandgenomics.rbind.io/post/understanding-p-value-multiple-comparisons-fdr-and-q-value/

- ROC Curves, https://www.r-bloggers.com/roc-curves/, https://github.com/dariyasydykova/open_projects/tree/master/ROC_animation

- The p value and the base rate fallacy, https://www.statisticsdonewrong.com/p-value.html

- Simulation-Based Power Analysis for ANOVA Designs, https://github.com/Lakens/ANOVApower

- `misc/week09.ppt` - Linear and Nonlinear Modeling, http://statwww.epfl.ch/davison/teaching/Microarrays/sched.html

- Ignatiadis, Nikolaos, Bernd Klaus, Judith B. Zaugg, and Wolfgang Huber. “Data-Driven Hypothesis Weighting Increases Detection Power in Genome-Scale Multiple Testing.” Nature Methods 13, no. 7 (2016): 577–80. https://doi.org/10.1038/nmeth.3885. - IHW - independent hypothesis weighting, a method for multiple testing correction weighting p-values by covariates, such as sum of read counts for gene expression (higher expression genes are prioritized), eQTLs (distance between variant and gene). Histograms of p-values, stratified by a covariate, indicates whether a covariate is informative. http://bioconductor.org/packages/release/bioc/html/IHW.html. Uses ideas from C. R. Genovese, K. Roeder, and L. Wasserman, “False Discovery Control with P-Value Weighting,” Biometrika 93, no. 3 (September 1, 2006): 509–24, https://doi.org/10.1093/biomet/93.3.509.