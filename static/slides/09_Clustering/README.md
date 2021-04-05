- `Clustering` - hierarchical clustering, distances, linkages
  - `lab`: `Average_linkage.R`, `Single_Linkage.R`, `Complete_Linkage.R`, `Divisive_Kmeans.R`, `CramersV_gower.R`
- `Clustering1` - K-means, SOM, PCA, NMF, t-SNE
  - `lab`: `Kmeans_PAM.R` (uses `BreastCancer.RData`), `PCA.R` (uses `nci60.tsv`), `MDS.Rmd`, 

## `hw` - homework material

## `lab` - in-class lab material

- `PCA.R` - uses `nci60.tsv` to demo PCA, batch, MDS, t-SNE, 3D PCA.
    - dist()
    - hclust()
    - distribution boxplots, filtering
    - PCA, with batch
    - MDS
    - t-SNE
    - UMAP
    - Differential expression with limma, heatmap
    - Noise - differential expression and clustering
    

- `NCI60Example` - PCA on NCI60 subset, accompanying Chen Meng et al., “Dimension Reduction Techniques for the Integrative Analysis of Multi-Omics Data,” Briefings in Bioinformatics 17, no. 4 (July 2016): 628–41, doi:10.1093/bib/bbv108. https://github.com/aedin/NCI60Example
- [nci60.tsv](http://odin.mdacc.tmc.edu/~kdo/TeachBioinf/Projects%20&%20Data%20Sets/nci60.tsv) - cell types can be clustered


## `read` - reading assignment

- Daxin Jiang, Chun Tang, and Aidong Zhang. “[Cluster Analysis for Gene Expression Data: A Survey](https://doi.org/10.1109/TKDE.2004.68).” IEEE Transactions on Knowledge and Data Engineering, (November 2004) - Clustering overview for gene expression studies. Definitions, proximity measures (Euclidean, Pearson), clustering (K-means, SOM, hierarchical, graph-theoretical, model-based, density, the use of PCA), biclustering. Metrics for clustering QC (homogeneity, separation, Rand, Jaccard, reliability)

- Patrik D’haeseleer, “[How Does Gene Expression Clustering Work?](https://doi.org/10.1038/nbt1205-1499),” Nature Biotechnology, (December 2005) - Clustering distances. Recommendations for gene expression choices of clustering

- Satagopan, Jaya M., and Katherine S. Panageas. “[A Statistical Perspective on Gene Expression Data Analysis](https://doi.org/10.1002/sim.1350).” Statistics in Medicine, (February 15, 2003) - Intro into microarray technology, statistical questions. Hierarchical clustering - clustering metrics. MDS algorithm. Class prediction - linear discriminant analysis algorithm and cross-validation. SAS and S examples

- Altman, Naomi, and Martin Krzywinski. “[Points of Significance: Clustering](https://doi.org/10.1038/nmeth.4299).” Nature Methods, (May 30, 2017) - Clustering depends on gene scaling, clustering method, number of simulations in k-means clustering.

- Krzywinski, Martin, and Naomi Altman. “[Points of Significance: Importance of Being Uncertain](https://doi.org/10.1038/nmeth.2613).” Nature Methods 10, no. 9 (September 2013)

- Altman, Naomi, and Martin Krzywinski. “[Points of Significance: Association, Correlation and Causation](https://doi.org/10.1038/nmeth.3587).” Nature Methods, (September 29, 2015)

## Dimensionality reduction

- Abdi, Hervé, and Lynne J. Williams. “[Principal Component Analysis](https://doi.org/10.1002/wics.101).” Wiley Interdisciplinary Reviews: Computational Statistics, (July 2010) - PCA in-depth review. Mathematical formulations, terminology, examples, interpretation. Figures showing PC axes, rotations, projections, circle of correlation. Rules for selecting number of components. Rotation - varimax, promax, illustrated. Correspondence analysis for nominal variables, Multiple Factor Analysis for a set of observations described by several groups (tables) of variables. Appendices - eigenvalues and eigenvectors, positive semidefinite matrices, SVD

- Wall, Michael. “[Singular Value Decomposition and Principal Component Analysis](https://link.springer.com/chapter/10.1007/0-306-47815-3_5),” - SVD and PCA statistical intro. Relation of SVD to PCA, Fourier transform. Examples of applications, including genomics. 

- Lever, Jake, Martin Krzywinski, and Naomi Altman. “[Points of Significance: Principal Component Analysis](https://doi.org/10.1038/nmeth.4346).” Nature Methods, (June 29, 2017) PCA explanation, the effect of scale. Limitations

- Lee, D. D., and H. S. Seung. “[Learning the Parts of Objects by Non-Negative Matrix Factorization](https://doi.org/10.1038/44565).” Nature, (October 21, 1999) - Non-negative matrix factorization (NMF) principles, compared with vector quantization (VQ) and PCA. Intuition behind NMF learning parts and PCA learning the whole.

- Lee, Daniel D., and H. Sebastian Seung. “[Algorithms for Non-Negative Matrix Factorization](http://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization.pdf).” In Advances in Neural Information Processing Systems, MIT Press, 2001 - Two algorithms for solving NMF - Euclidean distance and Kullback-Leibler divergence, with proofs.

- Meng, Chen, Oana A. Zeleznik, Gerhard G. Thallinger, Bernhard Kuster, Amin M. Gholami, and Aedín C. Culhane. “[Dimension Reduction Techniques for the Integrative Analysis of Multi-Omics Data](https://doi.org/10.1093/bib/bbv108).” Briefings in Bioinformatics, (July 2016) - Dimensionality reduction techniques - PCA and its derivatives, NMF. Table 1 - Terminology. Table 2 - methods, tools, visualization packages. Methods for integrative data analysis of multi-omics data.

- Lee, Su-In, and Serafim Batzoglou. “[Application of Independent Component Analysis to Microarrays](https://doi.org/10.1186/gb-2003-4-11-r76).” Genome Biology 4, no. 11 (2003) - Independent Components Analysis theory and applications.

- Stein-O’Brien, Genevieve L., Raman Arora, Aedin C. Culhane, Alexander V. Favorov, Lana X. Garmire, Casey S. Greene, Loyal A. Goff et al. "[Enter the matrix: factorization uncovers knowledge from omics](https://doi.org/10.1016/j.tig.2018.07.003)." Trends in Genetics, (2018) - Matrix factorization and visualization. Refs to various types of MF methods. Terminology, Fig 1 explanation of MF in terms of gene expression and biological processes. References to biological examples.

- Yeung, K. Y., and W. L. Ruzzo. “[Principal Component Analysis for Clustering Gene Expression Data](https://doi.org/10.1093/bioinformatics/17.9.763).” Bioinformatics, (September 1, 2001) - PCA is not always good for denoising data before clustering, clustering of PCs often worse than the original data. Simulated and real-life data. Data used for benchmarks: http://faculty.washington.edu/kayee/pca/

- Libbrecht, Maxwell W., and William Stafford Noble. “[Machine Learning Applications in Genetics and Genomics](https://doi.org/10.1038/nrg3920).” Nature Reviews. Genetics, (June 2015) - Machine learning in genomics. Supervised/unsupervised learning, semi-supervised, bayesian (incorporating prior knowledge), feature selection, imbalanced class sizes, missing data, networks. 

- Meng, Chen, Bernhard Kuster, Aedín C. Culhane, and Amin Moghaddas Gholami. “[A Multivariate Approach to the Integration of Multi-Omics Datasets](https://doi.org/10.1186/1471-2105-15-162).” BMC Bioinformatics (May 29, 2014) - MCIA - multiple correspondence analysis for integrating multiple datasets. Statistics and implementation in [omicade4](https://bioconductor.org/packages/omicade4/) - Multiple co-inertia analysis of omics datasets.

- [Singular Value Decomposition (SVD) Tutorial: Applications, Examples, Exercises](https://blog.statsbot.co/singular-value-decomposition-tutorial-52c695315254) blog post

- Liu, Yanchi, Zhongmou Li, Hui Xiong, Xuedong Gao, and Junjie Wu. “[Understanding of Internal Clustering Validation Measures](https://doi.org/10.1109/ICDM.2010.35),” IEEE, 2010 - Internal clustering validation metrics, table, concise description of each.

- Guido Kraemer, Markus Reichstein, and Miguel D. Mahecha, “[DimRed and CoRanking Unifying Dimensionality Reduction in R](https://doi.org/10.32614/RJ-2018-039),” The R Journal, 2018 - R packages implementing 15 methods for dimensionality reduction, from PCA, ICA, MDS to Laplasian eigenmaps. Brief but very good overview of each method, its complexity. Quality metrics to judge the quality of embedding. [GitHub](https://github.com/gdkrmr/dimRed)

- Amir, El-ad David, Kara L Davis, Michelle D Tadmor, Erin F Simonds, Jacob H Levine, Sean C Bendall, Daniel K Shenfeld, Smita Krishnaswamy, Garry P Nolan, and Dana Pe’er. “[ViSNE Enables Visualization of High Dimensional Single-Cell Data and Reveals Phenotypic Heterogeneity of Leukemia](https://doi.org/10.1038/nbt.2594).” Nature Biotechnology, (June 2013) - viSNE paper - tSNE (Barnes-Hut) implementation for single-cell data, and the `cyt` tool for visualization. [Supplementary methods](https://media.nature.com/original/nature-assets/nbt/journal/v31/n6/extref/nbt.2594-S1.pdf) - details of t-SNE algorithm, [Details of usage](https://www.denovosoftware.com/site/manual/visne.htm) 

- Belacel, Nabil, Qian Wang, and Miroslava Cuperlovic-Culf. “[Clustering Methods for Microarray Gene Expression Data](https://doi.org/10.1089/omi.2006.10.507).” Omics: A Journal of Integrative Biology, (2006) - Clustering methods overview. Hierarchical (agglomerative, divisive), partitional clustering (K-means, K-medoids, SOM). DBSCAN and other density-based algorithms. Graph-theoretical cllustering. Fuzzy clustering, expectation-maximization methods. Table with software.

- Kossenkov, Andrew V., and Michael F. Ochs. “[Matrix Factorisation Methods Applied in Microarray Data Analysis](https://doi.org/10.1504/IJDMB.2010.030968).” International Journal of Data Mining and Bioinformatics, (2010) - Matrix factorization methods for genomics data. SVD, PCA, ICA, NCA, NMF (sparse and least squares NMF), Bayesian decomposition

- Chavent, Marie, Vanessa Kuentz-Simonet, Amaury Labenne, and Jérôme Saracco. “[Multivariate Analysis of Mixed Data: The R Package PCAmixdata](http://arxiv.org/abs/1411.4911).” ArXiv, December 8, 2017 - PCAmixdata - R package for PCA on a mixture of numerical and categorical variables. Other packages - ade4, FactoMineR. Theory, statistics, code examples with interpretation. [PCAmixdata](https://CRAN.R-project.org/package=PCAmixdata)

### Tools

- [NbClust](https://CRAN.R-project.org/package=NbClust) - Determining the Best Number of Clusters in a Data Set. It provides 30 indexes for determining the optimal number of clusters in a data set and offers the best clustering scheme from different results to the user.

- [philentropy](https://CRAN.R-project.org/package=philentropy) - Similarity and Distance Quantification Between Probability Functions. Computes 46 optimized distance and similarity measures for comparing probability functions.

- [fpc](https://CRAN.R-project.org/package=fpc) - Flexible Procedures for Clustering R package. `prediction.strength` - function to calculate the optimal number of clusters.



## `misc` - misc presentations and materials

- `clustering.gif` - animation of agglomerative clustgering. [Source](https://cdn-images-1.medium.com/max/800/1*ET8kCcPpr893vNZFs8j4xg.gif)

- Clustering categorical variables using Jaccard and Tanimoto distances [https://stats.stackexchange.com/questions/49453/calculating-jaccard-or-other-association-coefficient-for-binary-data-using-matri](https://stats.stackexchange.com/questions/49453/calculating-jaccard-or-other-association-coefficient-for-binary-data-using-matri)

- Distances in `vegdist` function from `vegan` package [http://www.pmc.ucsc.edu/~mclapham/Rtips/cluster.htm](http://www.pmc.ucsc.edu/~mclapham/Rtips/cluster.htm)

- `n.dimensionreduction.pdf` - Dimension reduction. http://www.biostat.jhsph.edu/~iruczins/teaching/kogo/notes/n.dimensionreduction.pdf

- `lecture_08.pdf` - hierarchical clustering, agglomerative, steps. Ward's method. Single link, complete. [Source](http://www.stat.cmu.edu/~cshalizi/350/)

- `lecture-10.pdf` - Principal Components I. Reading assignment. http://www.stat.cmu.edu/~cshalizi/350/lectures/10/lecture-10.pdf

- `luxburg_clustering_lectures.pdf` - Lectures on Clustering by Ulrike von Luxburg, ~3.5 hours [http://videolectures.net/bootcamp07_luxburg_clu/#](http://videolectures.net/bootcamp07_luxburg_clu/#)

- [CRAN overview of clustering functions](http://cran.at.r-project.org/web/views/Cluster.html)

- PCA statistics explanation [1](http://users.ics.aalto.fi/jhollmen/dippa/node30.html), [2](https://onlinecourses.science.psu.edu/stat505/node/51)

- Google scholar links https://scholar.google.com/scholar?hl=en&q=cluster+analysis&btnG

- ConsensusClusterPlus - cluster count and membership, https://www.bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html

- sigclust: Statistical Significance of Clustering, https://cran.r-project.org/web/packages/sigclust/index.html

- https://datascienceplus.com/finding-optimal-number-of-clusters/

- `kohonen_som.pdf` - Self-Organizing Map, https://eric.univ-lyon2.fr/~ricco/cours/slides/en/kohonen_som.pdf


## ToDo - todo list

- Kullback-Leibler, Jensen-Shannon divergence, https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence

- Ester, M. “A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise.,” n.d. https://ocs.aaai.org/Papers/KDD/1996/KDD96-037.pdf. - DBSCAN paper, density-based spatial clustering for two-dimensional data. R tutorial at http://www.sthda.com/english/articles/30-advanced-clustering/105-dbscan-density-based-clustering-essentials/

- DBSCAN clustering tutorial, [https://medium.com/nearist-ai/dbscan-clustering-tutorial-dd6a9b637a4b](https://medium.com/nearist-ai/dbscan-clustering-tutorial-dd6a9b637a4b)


- MCL Markov Clustering algorithm, https://micans.org/mcl/, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC101833/

- Get some presentation material, and use as reading assignment. Abdi, Hervé, and Lynne J. Williams. “Principal Component Analysis.” Wiley Interdisciplinary Reviews: Computational Statistics 2, no. 4 (July 2010): 433–59. doi:10.1002/wics.101. - PCA in-depth review. Mathematical formulations, terminology, examples, interpretation. Figures showing PC axes, rotations, projections, circle of correlation. Rules for selecting number of components. Rotation - varimax, promax, illustrated. Correspondence analysis for nominal variables, Multiple Factor Analysis for a set of observations described by several groups (tables) of variables. Appendices - eigenvalues and eigenvectors, positive semidefinite matrices, SVD

- Padilha, Victor A., and Ricardo J. G. B. Campello. “A Systematic Comparative Evaluation of Biclustering Techniques.” BMC Bioinformatics 18, no. 1 (January 23, 2017): 55. https://doi.org/10.1186/s12859-017-1487-1. - Experimental evaluation of 17 biclustering methods, in simulated and real settings, 5 different experimental scenarios. Description of biclustering types and each specific method, links.

- Pontes, Beatriz, Raúl Giráldez, and Jesús S. Aguilar-Ruiz. “Biclustering on Expression Data: A Review.” Journal of Biomedical Informatics 57 (October 2015): 163–80. https://doi.org/10.1016/j.jbi.2015.06.028. - Detailed review of biclustering algorithms, methods.

- Prelić, Amela, Stefan Bleuler, Philip Zimmermann, Anja Wille, Peter Bühlmann, Wilhelm Gruissem, Lars Hennig, Lothar Thiele, and Eckart Zitzler. “A Systematic Comparison and Evaluation of Biclustering Methods for Gene Expression Data.” Bioinformatics (Oxford, England) 22, no. 9 (May 1, 2006): 1122–29. https://doi.org/10.1093/bioinformatics/btl060. - Review of the performance of five biclustering algorithms, and BiMax algorithm. Benchmarking on a real gene expression data, quality judged by functional enrichment. Description of the algorithms in the supplemental. http://people.ee.ethz.ch/~sop/bimax/



# `Legendre`, Misc functions

- `corPerm.R` - (P. Legendre): Three functions to test the Pearson correlation coefficient by permutation.

- `hclust.PL.R` - an extended version of hclust in R, by P. Legendre. In this function, the Ward1 algorithm is implemented by method= "ward.D" and the Ward2 algorithm by method="ward.D2". That function is available on the web page http://numericalecology.com, in the section on R-language functions. [Source](http://adn.biol.umontreal.ca/~numericalecology/Rcode/)

- `ward.pdf` - Ward’s Hierarchical Agglomerative Clustering Method: Which Algorithms Implement Ward’s Criterion? Fionn Murtagh , Pierre Legendre [PDF](https://link.springer.com/article/10.1007/s00357-014-9161-z)

- `MIT7_91JS14_Lecture8.pdf` - PCA https://www.youtube.com/watch?v=MniYgsZSp30&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=8

[Hierarchical, K-means, PAM, PCA demo](https://github.com/jennybc/stat540_2014/blob/master/seminars/seminar09_clustering-pca.rmd)



# Video

[Hierarchical clustering playlist](https://www.youtube.com/playlist?list=PLBv09BD7ez_7qIbBhyQDr-LAKWUeycZtx), and other clustering videos from [Viktor Lavrenko](https://www.youtube.com/user/victorlavrenko/playlists)

[Multidimensional Scaling, J.B. Kruskal. AT&T Bell Laboratories (1962)](http://stat-graphics.org/movies/multidimensional-scaling.html)

[Singular value decomposition](https://www.youtube.com/watch?v=YKmkAoIUxkU)

- t-SNE, https://www.youtube.com/watch?v=EMD106bB2vY

# ToDo

- [ClusterEnG](http://education.knoweng.org/clustereng/) - Interactive education in clustering

- https://github.com/gabrielodom/Bioc2019pathwayPCA - Workshop repository for pathwayPCA at the 2019 Bioconductor conference, https://gabrielodom.github.io/pathwayPCA/articles/

- Reese, Sarah E., Kellie J. Archer, Terry M. Therneau, Elizabeth J. Atkinson, Celine M. Vachon, Mariza de Andrade, Jean-Pierre A. Kocher, and Jeanette E. Eckel-Passow. “A New Statistic for Identifying Batch Effects in High-Throughput Genomic Data That Uses Guided Principal Component Analysis.” Bioinformatics (Oxford, England) 29, no. 22 (November 15, 2013): 2877–83. https://doi.org/10.1093/bioinformatics/btt480. - gPCA R package for quantifying the level of batch effect in the data. Introduction into PCA solved using SVD, guided PCA (input matrix as a product of batch and true matrix), test statistics quantifying the proportion of variance due to batch effect, its significance estimated using permutations.https://cran.r-project.org/web/packages/gPCA/index.html

- Clustering validity: "clValid", "NbClust". https://davetang.org/muse/2019/01/23/the-golden-rule-of-bioinformatics

- Some Technical Notes on Kullback-Leibler Divergence. https://alexanderetz.com/2018/09/05/some-technical-notes-on-kullback-leibler-divergence/

- Shlens, Jonathon. “A Tutorial on Principal Component Analysis.” ArXiv Preprint ArXiv:1404.1100, 2014. - PCA, SVD - detailed tutorial

- PCA explanation, math and R, http://www.hcbravo.org/IntroDataSci/lecture-note/pca/pca/, https://github.com/hcorrada/IntroDataSci/tree/master/materials/lectures/Unsupervised

Blashfield, Roger K., Mark S. Aldenderfer, and Leslie C. Morey. “11 Cluster Analysis Software.” In Handbook of Statistics, 2:245–66. Elsevier, 1982. https://doi.org/10.1016/S0169-7161(82)02014-8. - Blashfield et al report on a questionnaire regarding the use of cluster software, with fifty-three respondents yielding fifty different programs and packages. It is sometimes said that there are as many cluster methods as there are cluster analysis users.

Cluster Analysis is the mathematical study of methods for recognizing natural groups within a class of entities.

- Enright, A. J., S. Van Dongen, and C. A. Ouzounis. “An Efficient Algorithm for Large-Scale Detection of Protein Families.” Nucleic Acids Research 30, no. 7 (April 1, 2002): 1575–84. - Markov clustering to build networks. An adjacency matrix of -log10(p-values). https://micans.org/mcl/

- Yang, Jaewon, and Jure Leskovec. “Patterns of Temporal Variation in Online Media,” 177. ACM Press, 2011. https://doi.org/10.1145/1935826.1935863. - Time series clustering. The KSC algorithm that clusters time series with similar shapes using a metric that is invariant to scaling and shifting. http://s3l.stanford.edu/blog/?p=127

- List of software packages for multi-omics analysis. https://github.com/mikelove/awesome-multi-omics

- The indices comprise the Davies-Bouldin (DB) index [19], the Bayesian information criterion (BIC) [20], the silhouette index [21], the Calinski-Harabasz (CH) index [22], the Ball and Hall (BH) index [23], the Xu index [24], and the within-between (WB) index [25]. The DB index is a widely used and frequently cited whole-partition cluster validation index. (From Abu-Jamous and Kelly, “Clust.”)

- clustering scATACseq data: the TF-IDF way, https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/, Dimensionality Reduction for scATAC Data, http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/

- CCA aims to find linear combinations of features across data sets that are maximally correlated, identifying shared correlation structures across data sets. 
    - Butler, Andrew, Paul Hoffman, Peter Smibert, Efthymia Papalexi, and Rahul Satija. “Integrating Single-Cell Transcriptomic Data across Different Conditions, Technologies, and Species.” Nature Biotechnology 36, no. 5 (June 2018): 411–20. https://doi.org/10.1038/nbt.4096. - Meta-analysis of multiple scRNA-seq datasets, analogous to batch effect removal. Canonical Correlation Analysis (CCA). Learns a shared gene correlation structure that is conserved between the datasets using canonical correlation analysis. Works across technologies and organisms. Good methods description. Implemented in Seurat package.
    
- `CEMiTool` - gene co-expression analysis, reimplements WGCNA, includes selection of a soft-thresholding power using Cauchi distribution, gene enrichment analysis and, optionally, PPI network. Good overview of WGCNA algorithm. https://bioconductor.org/packages/release/bioc/html/CEMiTool.html
    - Russo, Pedro S. T., Gustavo R. Ferreira, Lucas E. Cardozo, Matheus C. Bürger, Raul Arias-Carrasco, Sandra R. Maruyama, Thiago D. C. Hirata, et al. “CEMiTool: A Bioconductor Package for Performing Comprehensive Modular Co-Expression Analyses.” BMC Bioinformatics 19, no. 1 (20 2018): 56. https://doi.org/10.1186/s12859-018-2053-1.

- PCA, Luis Serrano, by hand and visuals, https://youtu.be/g-Hb26agBFg

- Machine Learning — Singular Value Decomposition (SVD) & Principal Component Analysis (PCA). https://medium.com/@jonathan_hui/machine-learning-singular-value-decomposition-svd-principal-component-analysis-pca-1d45e885e491

- A collection of PCA methods. https://www.bioconductor.org/packages/release/bioc/html/pcaMethods.html

- How to cluster in High Dimensions. https://towardsdatascience.com/how-to-cluster-in-high-dimensions-4ef693bacc6

- No True Effects in High Dimensions. https://towardsdatascience.com/no-true-effects-in-high-dimensions-1f56360182cd

- Select Features for OMICs Integration, Univariate vs. Multivariate Feature Selection, by Nikolay Oskolkov. https://towardsdatascience.com/select-features-for-omics-integration-511390b7e7fd

- Workshop: Dimension reduction with R, Saskia Freytag, 19/07/2019. http://rpubs.com/Saskia/520216

The general procedure for traditional k-means clustering is as follows:
1. Randomly select k points in the dataset and assign the initial location of centroids to those points. 2. Calculate the squared Euclidean distance of all points within the dataset from each centroid.
3. Find each point’s closest centroid by this distance metric, and assign it to that centroid’s group.
4. Recalculate the location of each centroid as the mean of all its group members, assigned at step 3. 5. Repeat steps 2-4 several times.
Traditional k-means clustering, however, has several drawbacks. The results of k-means clustering are entirely dependent on the initial location of the centroids; these centroids simply roll toward their nearest local point density over the iterations described above (Arthur and Vassilvitskii, 2007). Because of this property, random selection of points to initialize the location of centroids can result in the initialization of centroids very close together within the dataset, and can also leave large areas of the dataset without a centroid initialized. This can lead to either a single population being labeled as two populations, or two populations being called a single population, respectively.
However, techniques have been developed to address this issue, with perhaps the most popular being k-means++ (Arthur and Vassilvitskii, 2007). K-means++ simply modifies how centroids are initially seeded. Rather than being seeded randomly, centroids are initialized sequentially as follows:
1. Thefirstcentroidisinitializeduniformlyrandomlytoapointinthedataset;inthiscase,eachpointhasanequalchanceofbeing selected.
2. For each point, the squared Euclidean distance is calculated in relation to its closest previously seeded centroid.
3. Eachadditionalcentroidistheninitializedrandomlyaswell,butusingthesedistancesasaweightedprobability;thisprobability
weighting increases the chances of choosing a point that is not near a previously selected centroid location.
4. Steps 2-3 are repeated until the user-defined number of centroids (k) has been initialized.
5. Traditional K-means clustering is then performed.
This process decreases the probability of centroids being seeded close to each other. However, this technique only takes into consideration the distance of a point from its nearest previously selected centroid, and not its distance from all previously selected centroids.

- Number of PCs in PCA to retain:
(i) The Scree-plot gives Eigenvalues versus number of principal components. The point of change (the elbow of the curvature) in the figures, which distinguishes the number of principal components, is the highest percentage to be retained.
(ii)Kaiser’s rule retains all components with Eigenvalues greater than one, and is a way of measuring the common variance of variables.
Franklin SB, Gibson DJ, Robertson PA, Pohlmann JT, Fralish JS: Parallel Analysis: a method for determining significant principal components. J Veg Sci 1995, 6:99–106.

- NMF does not have a unique solution: if we have a solution (W,H), then ($WD$,$D^{−1}H$) is another solution with any diagonal matrix D with positive diagonal elements.

- Stacklies, Wolfram, Henning Redestig, Matthias Scholz, Dirk Walther, and Joachim Selbig. “PcaMethods--a Bioconductor Package Providing PCA Methods for Incomplete Data.” Bioinformatics (Oxford, England) 23, no. 9 (May 1, 2007): 1164–67. https://doi.org/10.1093/bioinformatics/btm069. - PCA, various implementation. Probabilistic, Bayesian PCA. Imputation methods (SVDimpute, LLSimpute). https://www.bioconductor.org/packages/release/bioc/html/pcaMethods.html

- Rdimtools - Dimension Reduction and Estimation Methods. Many different methods in one interface. https://github.com/kyoustat/Rdimtools, https://cran.r-project.org/web/packages/Rdimtools/index.html

- Li, Yifeng, and Alioune Ngom. “The Non-Negative Matrix Factorization Toolbox for Biological Data Mining.” Source Code for Biology and Medicine 8, no. 1 (2013): 10. https://doi.org/10.1186/1751-0473-8-10. - The NMF MATLAB toolbox implementing most common NMF methods. Description and brief statistics of each method. Explanation how NMF is used for clustering, feature extraction. Contrasting PCA and NMF. https://sites.google.com/site/nmftool/

- Karim, Md Rezaul, Oya Beyan, Achille Zappa, Ivan G Costa, Dietrich Rebholz-Schuhmann, Michael Cochez, and Stefan Decker. “Deep Learning-Based Clustering Approaches for Bioinformatics.” Briefings in Bioinformatics, February 2, 2020, bbz170. https://doi.org/10.1093/bib/bbz170. - Deep learning for clustering. Overview of standard clustering approaches, their limitations. Benchmarked on clustering of images (BACH, breast cancer histology images https://iciar2018-challenge.grand-challenge.org/Home/), drug reviews and 10-star rating https://www.drugs.com/, several TCGA gene expression datasets. Convolutional autoencoders and density-based clustering perform well. GitHub repo with notebooks implementing benchmarks and links to papers and repos of the deep learning clustering methods https://github.com/rezacsedu/Deep-learning-for-clustering-in-bioinformatics

- 10 Clustering Algorithms With Python, https://machinelearningmastery.com/clustering-algorithms-with-python/

- Cluster Validation Statistics: Must Know Methods, https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/

- MDS: we use the classical MDS which finds the p dimensional embedding
X = [x1, x2, · · · , xn]T ∈ Rn×p that minimizes the loss function: loss = ∥XXT − G∥F , where
G=−1(In−111T)D(In−111T),Dij isthedistancebetweeni-thandj-thcells,and1denotes 2nn
a column vector of all ones. - from Li, Xinjun, Fan Feng, Wai Yan Leung, and Jie Liu. “ScHiCTools: A Computational Toolbox for Analyzing Single-Cell Hi-C Data.” Preprint. Bioinformatics, September 18, 2019. https://doi.org/10.1101/769513.
