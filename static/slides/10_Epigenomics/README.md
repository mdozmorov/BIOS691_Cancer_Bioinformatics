## `hw` - homework material

## `lab` - in-class lab material

- `Enrichment_Regions.Rmd` - Peak annotation and functional enrichment

## References

### Epigenomics intro

- Stephen B. Baylin and Peter A. Jones, “[A Decade of Exploring the Cancer Epigenome — Biological and Translational Implications](https://doi.org/10.1038/nrc3130),” Nature Reviews Cancer, (September 23, 2011) - Cancer epigenomics introduction, therapies

- Kagohara, Luciane T., Genevieve L. Stein-O’Brien, Dylan Kelley, Emily Flam, Heather C. Wick, Ludmila V. Danilova, Hariharan Easwaran, et al. “[Epigenetic Regulation of Gene Expression in Cancer: Techniques, Resources and Analysis](https://doi.org/10.1093/bfgp/elx018).” Briefings in Functional Genomics, August 11, 2017 - Review of epigeneitc modifications, methylation, histones, chromatin states, 3D. Technologies, databases, software. Lots of references

- Li, Bing, Michael Carey, and Jerry L. Workman. “[The Role of Chromatin during Transcription](http://www.cell.com/abstract/S0092-8674(07)00109-2).” Cell, (February 23, 2007) - Transcription process and the role of chromatin modifications. 

- Zhou, Vicky W., Alon Goren, and Bradley E. Bernstein. “[Charting Histone Modifications and the Functional Organization of Mammalian Genomes](https://www.nature.com/articles/nrg2905).” Nature Reviews. Genetics, (January 2011) - Histone marks review, ChIP-seq. Graphics of histone marks roles

- Wang, Zhibin, Dustin E. Schones, and Keji Zhao. “[Characterization of Human Epigenomes](https://doi.org/10.1016/j.gde.2009.02.001).” Current Opinion in Genetics & Development, (April 2009) - Concise description of main histone marks, their roles in transcription, and the corresponding studies. Figure 2 - schematic distribution of histone marks with respect to genes-TSSs.

- Zhang, Z. D., A. Paccanaro, Y. Fu, S. Weissman, Z. Weng, J. Chang, M. Snyder, and M. B. Gerstein. “[Statistical Analysis of the Genomic Distribution and Correlation of Regulatory Elements in the ENCODE Regions](https://doi.org/10.1101/gr.5573107).” Genome Research, (June 1, 2007) - ENCODE pilot project analysis. Non-random location of regulatory elements. Enrichment in TSSs, not in the middle or end of transcription sites. PCA and biplot representation of interrelatedness among TFs and histone marks, clustering.

- Huen, David S., and Steven Russell. “[On the Use of Resampling Tests for Evaluating Statistical Significance of Binding-Site Co-Occurrence](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-359).” BMC Bioinformatics (June 30, 2010) - Thorough review on permutation. Problem with different numbers of ROIs/epiregions in permutations. Regions may overlap (i.e., occur in clusters). Independent assignment is best.

### Epigenomic enrichment

- McLean, Cory Y., Dave Bristor, Michael Hiller, Shoa L. Clarke, Bruce T. Schaar, Craig B. Lowe, Aaron M. Wenger, and Gill Bejerano. “[GREAT Improves Functional Interpretation of Cis-Regulatory Regions](https://doi.org/10.1038/nbt.1630).” Nature Biotechnology 28, no. 5 (May 2010) - Hypergeometric and binomial enrichment of regulatory regions in relation to genesgenomic regions and their ontologies.

- Dozmorov, Mikhail G. “[Epigenomic Annotation-Based Interpretation of Genomic Data: From Enrichment Analysis to Machine Learning](https://doi.org/10.1093/bioinformatics/btx414).” Bioinformatics (Oxford, England) 33, no. 20 (October 15, 2017)

### Hidden Markov Models, genome segmentation overview

- Eddy, Sean R. “[What Is a Hidden Markov Model?](https://doi.org/10.1038/nbt1004-1315)” Nature Biotechnology, (October 2004)

- Eddy, S. R. “[Multiple Alignment Using Hidden Markov Models](https://pdfs.semanticscholar.org/bc1e/ffe17026c1396697e4ea2b399f7a049202fc.pdf).” Proceedings. International Conference on Intelligent Systems for Molecular Biology (1995) - HMMER hidden markov models for protein sequence alignment. First paper

- Schuster-Böckler, Benjamin, and Alex Bateman. “[An Introduction to Hidden Markov Models](https://doi.org/10.1002/0471250953.bia03as18).” 2007 - Hidden Markov Models primer

- L.R. Rabiner, “[A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition](https://doi.org/10.1109/5.18626),” Proceedings of the IEEE, (February 1989) - Theory and statistics of Hidden Markov Models. Very detailed, thorough

- Ernst, Jason, and Manolis Kellis. “[Chromatin-State Discovery and Genome Annotation with ChromHMM](https://doi.org/10.1038/nprot.2017.124).” Nature Protocols, (December 2017) - ChromHMM protocol. Intro about ChromHMM, other methods. Links to genome annotation models.

- Hoffman, Michael M., Orion J. Buske, Jie Wang, Zhiping Weng, Jeff A. Bilmes, and William Stafford Noble. “[Unsupervised Pattern Discovery in Human Chromatin Structure through Genomic Segmentation](https://doi.org/10.1038/nmeth.1937).” Nature Methods 9, no. 5 (May 2012) - Segway - segmentation and prediction of genomic states. Using dynamic Bayesian network

<!--
- Rojano, Elena, Pedro Seoane, Juan A G Ranea, and James R Perkins. “Regulatory Variants: From Detection to Predicting Impact.” Briefings in Bioinformatics, June 8, 2018. [https://doi.org/10.1093/bib/bby039](https://doi.org/10.1093/bib/bby039) - Regulatory SNPs methods of analysis review. Key regulatory elements (promoters, enhancers, silencers, insulators). Clustering and redundancy of regulatory elements. Post-transcriptional regulation (miRNA, lincRNA). How to call SNPs. Table 1 - databases of regulatory elements. Annotations for predicting pathogenicity (CADD, DANN, GWAVA, FATHMM-MKL, LINSIGHT, others, Table 3). Table 2 - annotation tools. Table 4 - current experimental techniques, including Hi-C. Example diseases having regulatory variants.
-->

## `misc` - misc presentations and materials

- `bw_doodling_markov.pdf` - A simple Markov model in R

- `hidden_markov_models.pdf` - Hidden Markov Models. Ben Langmead. http://www.cs.jhu.edu/~langmea/resources/lecture_notes/hidden_markov_models.pdf

- `lecture_17.pdf` - Hidden Markov Model reading. https://ocw.mit.edu/courses/mathematics/18-417-introduction-to-computational-molecular-biology-fall-2004/lecture-notes/lecture_17.pdf

- `markov_chains.pdf` - Markov chain for CpG discovery. Ben Langmead, 1st lec http://www.cs.jhu.edu/~langmea/resources/lecture_notes/markov_chains.pdf

- `MIT7_91JS14_Lecture10.pdf` - Markov & Hidden Markov Models of Genomic & Protein Features. https://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/lecture-slides/MIT7_91JS14_Lecture10.pdf, https://www.youtube.com/watch?v=d5NMrA2HkG4

## ToDo - todo list

- Segway dynamic Bayesian network explanation - ~15min https://www.youtube.com/watch?v=iKLvCuFD1MA&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=18, and at Hoffman, Michael M., Orion J. Buske, Jie Wang, Zhiping Weng, Jeff A. Bilmes, and William Stafford Noble. “Unsupervised Pattern Discovery in Human Chromatin Structure through Genomic Segmentation.” Nature Methods 9, no. 5 (May 2012): 473–76. doi:10.1038/nmeth.1937.


- http://setosa.io/ev/markov-chains/

- A friendly introduction to Bayes Theorem and Hidden Markov Models. https://www.youtube.com/watch?v=kqSzLo9fenk&list=UUgBncpylJ1kiVaPyP-PZauQ&index=17, and the Python notebook, https://github.com/luisguiserrano/hmm

- IDEAS chromatin segmentation for 127 cell lines from Roadmap Epigenomics. Critique of ChromHMM, proposing Bayesian non-parametric technique to automatically choose the number of states from the data. Uses signal, not binarized. Data sources, include promoter-capture HiC data. http://personal.psu.edu/yzz2/IDEAS/, https://main.genome-browser.bx.psu.edu/
    - Zhang, Yu, and Ross C. Hardison. “Accurate and Reproducible Functional Maps in 127 Human Cell Types via 2D Genome Segmentation.” Nucleic Acids Research 45, no. 17 (September 29, 2017): 9823–36. https://doi.org/10.1093/nar/gkx659.

