+ `01_Alignment_intro.Rmd` - Alignment introduction

+ `lab/naive_exact.R` - naive exact matching  
- `misc/seq_aln.pdf` - Sequencing technology, (from slide 37) local alignment, Smith-Waterman algorithm. http://www.biostat.jhsph.edu/~khansen/teaching/2014/140.668/seq_aln.pdf
+ `02_Approximate_matching_Langmead.pdf` - Hamming/Edit distance, global/local alignment overview
+ `lab/naive_Hamming.R` - Hamming distance matching  
+ `lab/edDistRecursive.R` - edit distance, recursive
+ `lab/edDistDynamic.R` - edit distance, dynamic programming  


- `03_Needleman_Wunsch.Rmd` - Needleman-Wunsch global alignment basics
- `04_Smith-Waterman.pdf` - Smith-Waterman algorithm, after Ben Langmead http://www.biostat.jhsph.edu/~khansen/teaching/2014/140.668/seq_aln.pdf
- `05_BWT.pdf` - Short overview of BWT, from MIT7_91JS14_Lecture5.pdf, https://www.youtube.com/watch?v=P3ORBMon8aw&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=5

- `17 - Genome arithmetic with bedtools.pdf` - BEDtools tutorial. From https://github.com/quinlan-lab/applied-computational-genomics
- `18 - Real world analyses with bedtools..pdf` - BEDtools exercise, overlap of GWAS snps with enhancers, calculating significance of overlaps. From https://github.com/quinlan-lab/applied-computational-genomics


## `lab` - in-class lab material

- [TeachEnG](http://teacheng.illinois.edu/#nw) - educational tool for aslgnment algorithms, plylogenetic tree.  - Kim, Minji, Yeonsung Kim, Lei Qian, and Jun S. Song. “[TeachEnG: A Teaching Engine for Genomics](https://doi.org/10.1093/bioinformatics/btx447).” Bioinformatics, October 15, 2017

- [kcmbt_mt](https://github.com/abdullah009/kcmbt_mt) - KCMBT: A very fast k-mer counter

- [KMC](http://sun.aei.polsl.pl/kmc. https://github.com/refresh-bio/KMC) - counting k-mers from (possibly gzipped) FASTQ/FASTA files


## References

- Nagarajan, Niranjan, and Mihai Pop. “[Sequence Assembly Demystified](https://doi.org/10.1038/nrg3367).” Nature Reviews. Genetics, March 2013 - Gentle introduction into genome assembly. Technologies. Box2: Greedy, overlap-layout-consensus, De Bruijn. Problems

- Pevzner, P. A., H. Tang, and M. S. Waterman. “[An Eulerian Path Approach to DNA Fragment Assembly](https://doi.org/10.1073/pnas.171285098).” Proceedings of the National Academy of Sciences, August 14, 2001 - First de Bruijn graph for genome assembly paper. Idea of breaking reads into fragments. Typical approach reads are vertices connected by edges if they overlap. Hamiltonian path problem - visit each vertex exactly once, NP-complete. de Bruijn graph - overlapping fragments are edges, and the problem is Eulerian path - visit each edge once. Error-correction algorithm.

- Compeau, Phillip E C, Pavel A Pevzner, and Glenn Tesler. “[How to Apply de Bruijn Graphs to Genome Assembly](https://doi.org/10.1038/nbt.2023).” Nature Biotechnology, November 8, 2011

- Chaisson, Mark J. P., Richard K. Wilson, and Evan E. Eichler. “[Genetic Variation and the de Novo Assembly of Human Genomes](https://doi.org/10.1038/nrg3933).” Nature Reviews Genetics, October 7, 2015 - Genome assembling strategies, problems. OLC, De Bruijn, string graphs. Types of gaps. 

- Miller, Jason R., Sergey Koren, and Granger Sutton. “[Assembly Algorithms for Next-Generation Sequencing Data](https://doi.org/10.1016/j.ygeno.2010.03.001).” Genomics, June 2010 - Assembly tools for overlap/layout/consensus and the de Bruijn graph approaches. de Bruin graph Issues with genome assembly, potential solutions.

- String Graph Assembler. Simpson, J. T., and R. Durbin. “[Efficient de Novo Assembly of Large Genomes Using Compressed Data Structures](http://www.genome.org/cgi/doi/10.1101/gr.126953.111).” Genome Research, March 1, 2012 - SGA - String Graph Assembler. From an FM-index. Velvet, ABySS, SOAPdenovo de Bruijn graph assemblers. BWA and FM explanation

- Koren, Sergey, and Adam M. Phillippy. “[One Chromosome, One Contig: Complete Microbial Genomes from Long-Read Sequencing and Assembly](https://doi.org/10.1016/j.mib.2014.11.014).” Current Opinion in Microbiology, February 2015 - Genome assembly overview focusing on long reads. Repeats (global and local) are problematic. Details on technologies: PacBio RS, Illumina's Moleculo, ONT MinION. Assembling approaches: OLC, hierarchical hybrid (long reads correction using another technology) and non-hybrid (self long reads alignment-correction). Assembly augmentation: gap filling, scaffolding, read threading. Table 1 - long read assembly tools and descriptions.

### Alignment algorithms and tools

- [textdistance](https://github.com/life4/textdistance) - compute distance between sequences. 30+ algorithms, pure python implementation, common interface, optional external libs usage

- Needleman, S. B., and C. D. Wunsch. “[A General Method Applicable to the Search for Similarities in the Amino Acid Sequence of Two Proteins](https://doi.org/10.1016/0022-2836(70)90057-4).” Journal of Molecular Biology, March 1970
    - [Online demo](http://experiments.mostafa.io/public/needleman-wunsch/)
    - [Wikipedia](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)

- Smith, T. F., and M. S. Waterman. “[Identification of Common Molecular Subsequences](https://doi.org/10.1016/0022-2836(81)90087-5).” Journal of Molecular Biology, March 25, 1981

- Burrows, Michael, and David J Wheeler. “[A Block-Sorting Lossless Data Compression Algorithm](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.8069),” 1994. - BWT paper

- Ferragina, Paolo, and Giovanni Manzini. “[Opportunistic Data Structures with Applications](https://dl.acm.org/citation.cfm?id=796543).” In Foundations of Computer Science, IEEE, 2000 - FM index paper

- Li, Heng, and Richard Durbin. “[Fast and Accurate Short Read Alignment with Burrows-Wheeler Transform](https://doi.org/10.1093/bioinformatics/btp324).” Bioinformatics, July 15, 2009 - BWA first paper

- Li, Heng, and Richard Durbin. “[Fast and Accurate Long-Read Alignment with Burrows-Wheeler Transform](https://doi.org/10.1093/bioinformatics/btp698).” Bioinformatics, March 1, 2010 - BWA second paper

- Li, Heng. “[Aligning Sequence Reads, Clone Sequences and Assembly Contigs with BWA-MEM](https://arxiv.org/abs/1303.3997).” ArXiv, 2013. - BWA-MEM (maximal exact matches). Automatically select end-to-end or gapped alignment.

- Langmead, Ben, and Steven L Salzberg. “[Fast Gapped-Read Alignment with Bowtie 2](https://doi.org/10.1038/nmeth.1923).” Nature Methods, March 4, 2012 - Bowtie2 gapped alignment

- Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. “[STAR: Ultrafast Universal RNA-Seq Aligner](https://doi.org/10.1093/bioinformatics/bts635).” Bioinformatics, January 1, 2013 - STAR gapped aligner. Algorithm, testing on simulated dataset.

- Kim, Daehwan, Ben Langmead, and Steven L. Salzberg. “[HISAT: A Fast Spliced Aligner with Low Memory Requirements](https://doi.org/10.1038/nmeth.3317).” Nature Methods, April 2015 - HISAT aligner. Two indexes - one large with many small

- Fonseca, Nuno A., Johan Rung, Alvis Brazma, and John C. Marioni. “[Tools for Mapping High-Throughput Sequencing Data](https://doi.org/10.1093/bioinformatics/bts605).” Bioinformatics, December 15, 2012 - Sequencing mappers, lists and characteristics of DNA, RNA, bisulfite aligners.

- [Bedtools Cheatsheet](https://gist.github.com/ilevantis/6d6ecf8718a5803acff736c2dffc933e#subtract)



### Videos

- [Algorithms for DNA Sequencing](https://www.youtube.com/playlist?list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA) - 55 short videos by Ben Langmead
    - [Ben Langmead's Slides for Algorithms for DNA Sequencing Coursera class, GitHub repo](https://github.com/BenLangmead/ads1-slides.git)
    - [Ben Langmead's teaching material with links to videos, lecture notes, python notebooks](http://www.langmead-lab.org/teaching-materials/)

- [MIT 7.91J Foundations of Computational and Systems](https://www.youtube.com/watch?v=lJzybEXmIj0&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=1) - full lecture videos (approx 1h 20m), by Christopher Burge, David Gifford, Ernest Fraenkel

- Vector illustrations of BWT based searching on small strings, https://github.com/ekg/drawbwt
    - [FM-Indexes and Backwards Search](https://alexbowe.com/fm-index/)


## `misc`

- `02_Approximate_matching_Langmead_full.pdf` - Hamming/Edit distance, global/local alignment overview, extended version

- `MIT7_91JS14_Lecture6.pdf` - Gifford. overlap/consensus/layout, shortest common superstring, greedy algorithm, de Bruijn graph. https://www.youtube.com/watch?v=ZYW2AeDE6wU&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=6. For `01_Alignment_intro.Rmd`

- `MIT7_91JS14_Lecture5.pdf` - Gifford. Negative Binomial, BWT, FM index, https://www.youtube.com/watch?v=P3ORBMon8aw&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=5. For `01_Alignment_intro.Rmd`

- `MIT7_91JS14_Lecture3.pdf` - Gifford. Global/local alignment, NW/SW algorithms, dynamic programming, from video https://www.youtube.com/watch?v=PdyARRNwi7I&list=PLUl4u3cNGP63uK-oWiLgO7LLJV6ZCWXac&index=3. For `01_Alignment_intro.Rmd`

- `02_Approximate_matching_Langmead.pdf`
    - "[Approximate matching, Hamming and edit distance](https://www.youtube.com/watch?v=MCvHeW13DsE&index=30&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA)", [naive_Hamming.R](../../slides/04_Alignment/lab/naive_Hamming.R)
    - "[Solving the edit distance problem](https://www.youtube.com/watch?v=8Q2IEIY2pDU&index=33&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA)", [edDistRecursive.R](../../slides/04_Alignment/lab/edDistRecursive.R)
    - "[Using dynamic programming for edit distance](https://www.youtube.com/watch?v=0KzWq118UNI&index=34&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA)", [edDistDynamic.R](../../slides/04_Alignment/lab/edDistDynamic.R)
    - "[Edit distance for approximate matching](https://www.youtube.com/watch?v=NjfNZzJiu8o&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA&index=36)"
    - "[Meet the family: global and local alignment](https://www.youtube.com/watch?v=-bjSPP2v6_Q&index=37&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA)"

- `assembly_dbg.pdf` - To use. Langmead. de Brujin graph lecture, https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf

- `Basics_of_alignment_lecture1.pdf` - Skiena. Excerpt from http://www.algorist.com/computational_biology/pdf/lecture1.pdf, video https://www.youtube.com/watch?v=pJAIalWDgY0&list=PLA48145CC64FE7990

- `lecture4.pdf` - Skiena. Mapping, suffix tries, trees, arrays, searching. http://www.algorist.com/computational_biology/pdf/lecture4.pdf

- `seq_aln.pdf` - Sequencing technology, alignment, Spith-Waterman algorithm. Kasper Hansen, http://www.biostat.jhsph.edu/~khansen/teaching/2014/140.668/seq_aln.pdf


## ToDo

- [alignment-and-variant-calling-tutorial](https://github.com/ekg/alignment-and-variant-calling-tutorial) - basic walk-throughs for alignment and variant calling from NGS sequencing data, by Erik Garrison. 

