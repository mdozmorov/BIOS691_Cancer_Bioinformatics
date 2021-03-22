- `08_functional_enrichment` - Functional enrichment slides
- `lab/Enrichment_Genes.Rmd` - differential expression using limma, annotations, functional enrichment
- `lab/EnrichmentOverlay.Rmd` - overlaying DEGs over a pathway


## `lab` - in-class lab material

In-class enrichment analysis quiz https://docs.google.com/forms/d/e/1FAIpQLScqOW3bmOij-FQacKSoCUEqTKcPcJz5P-4gTOkrkVA1OcKVjw/viewform?usp=sf_link#start=openform


- `gene_set_testing.Rmd` - another demo of differential expression and functional enrichment using mroast

## `read` - reading assignment

- Reimand, Jüri, Ruth Isserlin, Veronique Voisin, Mike Kucera, Christian Tannus-Lopes, Asha Rostamianfar, Lina Wadi, et al. “Pathway Enrichment Analysis and Visualization of Omics Data Using g:Profiler, GSEA, Cytoscape and EnrichmentMap.” Nature Protocols, January 21, 2019. https://doi.org/10.1038/s41596-018-0103-9. - Pathway enrichment analysis tutorial. Main tools - g:Profiler, GSEA, Cytoscape, EnrichmentMap. Definitions, step-by-step protocol.

- Bard, Jonathan B. L., and Seung Y. Rhee. “Ontologies in Biology: Design, Applications and Future Challenges.” Nature Reviews. Genetics 5, no. 3 (March 2004): 213–22. doi:10.1038/nrg1295. - Ontologies review

- Ackermann, Marit, and Korbinian Strimmer. “A General Modular Framework for Gene Set Enrichment Analysis.” BMC Bioinformatics 10, no. 1 (2009): 47. doi:10.1186/1471-2105-10-47. - All steps for enrichment analysis, methods, statistics, GSEA.

- Efron, Bradley, and Robert Tibshirani. “On Testing the Significance of Sets of Genes.” The Annals of Applied Statistics, 2007, 107–29. - maxmean statistics for enrichment analysis. Comparison with GSEA.

- Huang, Da Wei, Brad T. Sherman, and Richard A. Lempicki. “Bioinformatics Enrichment Tools: Paths toward the Comprehensive Functional Analysis of Large Gene Lists.” Nucleic Acids Research 37, no. 1 (January 2009): 1–13. doi:10.1093/nar/gkn923. - Gene enrichment analyses tools. Statistics, concept of background. 68 tools, table

- Hung, Jui-Hung, Tun-Hsiang Yang, Zhenjun Hu, Zhiping Weng, and Charles DeLisi. “Gene Set Enrichment Analysis: Performance Evaluation and Usage Guidelines.” Briefings in Bioinformatics 13, no. 3 (May 2012): 281–91. doi:10.1093/bib/bbr049. - Details of GSEA. Statistics, correction for multiple testing. Lack of gold standard - concept of mutual coverage.

- Khatri, Purvesh, Marina Sirota, and Atul J. Butte. “Ten Years of Pathway Analysis: Current Approaches and Outstanding Challenges.” PLoS Computational Biology 8, no. 2 (2012): e1002375. doi:10.1371/journal.pcbi.1002375. - Review of enrichment analyses techniques, focus on pathways, Table 1 lists tools, limitations.

- Hänzelmann, Sonja, Robert Castelo, and Justin Guinney. “GSVA: Gene Set Variation Analysis for Microarray and RNA-Seq Data.” BMC Bioinformatics 14 (January 16, 2013): 7. doi:10.1186/1471-2105-14-7. - GSVA - a GSE method that estimates variation of pathway activity over a sample population in an unsupervised manner

- Khatri, Purvesh, Marina Sirota, and Atul J. Butte. “Ten Years of Pathway Analysis: Current Approaches and Outstanding Challenges.” Edited by Christos A. Ouzounis. PLoS Computational Biology 8, no. 2 (February 23, 2012): e1002375. https://doi.org/10.1371/journal.pcbi.1002375. - Overview of pathway analysis methods (ORA - over-representation analysis, FCS - functional class scoring, PT - pathway topology) and tools. Methods, limitations.

- Dave’s blog (http://davetang.org/muse/) search for “Gene ontology enrichment analysis” 
- Nam D., and Seon-Young K.. “**Gene-Set Approach for Expression Pattern Analysis.**” _Briefings in Bioinformatics_ 2008 https://www.ncbi.nlm.nih.gov/pubmed/18202032
- Mutation Consequences and Pathway Analysis working group. “**Pathway and Network Analysis of Cancer Genomes.**” _Nature Methods_ 2015 https://www.ncbi.nlm.nih.gov/pubmed/26125594
- de Leeuw, C. et.al. “**The Statistical Properties of Gene-Set Analysis.**” _Nature Reviews_ 2016 https://www.ncbi.nlm.nih.gov/pubmed/27070863

## `misc` - misc presentations and materials

- `lecture21.pdf`, `code21.*` - Hypergeometric distribution, Fisher's exact test. http://www.biostat.jhsph.edu/~iruczins/teaching/140.652/lecture21.pdf

- `Lesson12.pdf` - stats for chi-square distribution, test, Fisher's exact, McNemar test, http://www.biostat.umn.edu/~susant/Fall09ph6414/Lesson12.pdf

- `Gene_list_enrichment_Mar10.pdf` - full gene set enrichment analysis lecture, http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf


## `ToDo`

- Fast Gene Set Enrichment Analysis, https://github.com/ctlab/fgsea

- DESeq results to fgsea, https://github.com/stephenturner/deseq-to-fgsea

- Afsari, Bahman, Donald Geman, and Elana J. Fertig. “Learning Dysregulated Pathways in Cancers from Differential Variability Analysis.” Cancer Informatics 13, no. Suppl 5 (2014): 61–67. https://doi.org/10.4137/CIN.S14066. - Detection of dysregulated pathways - over-representation, enrichment, differential variability. Focus on the latter. EVA - expression variation analysis, DIRAC. Intuitive statistical explanation. GSReg package.

- Geistlinger, Ludwig, Gergely Csaba, and Ralf Zimmer. “Bioconductor’s EnrichmentBrowser: Seamless Navigation through Combined Results of Set- & Network-Based Enrichment Analysis.” BMC Bioinformatics 17, no. 1 (December 2016): 45. https://doi.org/10.1186/s12859-016-0884-1. - EnrichmentBrowser - gene set- and network-based enrichment analysis and visualization. Input - text-based data into summarized experiment, normalization, differential expression, one or multiple enrichment analyses, visualization. https://bioconductor.org/packages/release/bioc/html/EnrichmentBrowser.html

- `PaintOmics 3` - web tool for KEGG pathway enrichment analysis and visualization of gene expression (also, metabolite, protein, region-based data) over pathway diagrams. Competitors: MapMan, KaPPA-View, Pathview Web. Auto-detection of IDs. Analyzes fold change, time course. http://www.paintomics.org/
    - Hernández-de-Diego, Rafael, Sonia Tarazona, Carlos Martínez-Mira, Leandro Balzano-Nogueira, Pedro Furió-Tarí, Georgios J. Pappas, and Ana Conesa. “PaintOmics 3: A Web Resource for the Pathway Analysis and Visualization of Multi-Omics Data.” Nucleic Acids Research 46, no. W1 (July 2, 2018): W503–9. https://doi.org/10.1093/nar/gky466.
