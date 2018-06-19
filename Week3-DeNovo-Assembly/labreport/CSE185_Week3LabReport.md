## Abstract  

In past labs, we aligned reads to a known reference genome. This week, we were given two paired-end libraries corresponding to DNA fragments of the bacteria *S. aureus*, and were expected to assemble a genome from these libraries. We were able to demonstrate a de novo assembly pipeline by completing kmer-based error correction of the long-frag libraries, contig assembly, scaffolding, and gap closing, while consistently evaluating the assembly quality and working towards a finished genome. Ultimately, while the genome aligned well to the reference genome, there were many local misassemblies, and we outlined steps to take in the future to improve assembly quality.  

## Introduction  

In [Nature Reviews](https://www.nature.com/articles/nrg3367), Nagarajan and Pop explore recent advancements in genome assembly, including the development of de novo assembly strategies and the effective use of complementary information such as mate-pair or paired-end data. In this lab, we leverage both advancements toward the creation of our bacterial genome. De novo assembly is the process of creating a genome without the aid of a reference genome. This process is used when the genome of the organism of interest in unknown, such as for non-model, divergent, or rare species. De novo assembly can also help biologists study large structural variations from the reference genome. Even without a complete genome, partial assembly can help in gene finding and comparative analyses across species.  

This week, we assembled the reads of bacteria *S. aureus* into a complete genome without the aid of the reference. We are using this organism because of the small genome and the known reference to which we can evaluate our assembly at the end of the report. We were given four raw data files corresponding to two paired-end libraries. The first has fragments of 180 bp with short insert size (101 bp), and the second has fragments of 37 bp with long insert size (3500 bp). The first library is used to assemble the contigs, and the second library is used to scaffold.  

## Methods  

We first inspect the raw data by checking the number of reads for each file and running  **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v.0.11.7)** on each file. Evidently, the data is pretty bad, so we move forward with a kmer-based error correction method for the short insert size library.  

We use **[Jellyfish](http://www.genome.umd.edu/jellyfish.html) (v.2.2.7)** to count number of times a 31-mer occurs, and number of 31-mers for each count for this library. Using the histo command from jellyfish as well as features in iPython, we can make a histogram of the frequencies. The plot shows a high frequency of 31-mers occurring only once, indicating sequencing errors, so we run error corrector commands to combat this. **[SOAPdenovo2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3626529/)** error modules **KmerFreq_HA (v.2.03)** and **Corrector_HA (v.2.03)** find high-frequency k-mers similar to the low-frequency k-mers and corrects the sequencing errors so we can still use the untrustworthy kmers in the assembly. We then output corrected FASTQ files for the short insert size files, and when we rerun Jellyfish, we can see that the sequencing errors are omitted and the histogram follows the expected Poisson distribution.  

With these new corrected files, we rerun some of our previous tests — we find the number of reads, valley and peak, average read length, total number of base pairs, and estimated genome size. Now, we can move onto assembly using **[Minia](http://minia.genouest.org/) (v.2.0.7)** with kmer size 31 and determine basic statistics representing the goodness of this assembly using the contigs file and the webtool **[QUAST](http://quast.bioinf.spbau.ru/)**.  

Using class data, we assess kmer size versus maximum contig length and N50 score to determine how kmer length choice influences success of assembly.  

Now, we can move onto using the short jump library to scaffold the contigs we have using **[SSPACE](https://github.com/nsoranzo/sspace_basic) (v.2.0)**. We assess the goodnes of the output by recording the number of scaffolds, the max scaffold length, and the N50. However, we know these files have many low-quality base calls from the FastQC analysis at the beginning of this lab. So, we rerun SSPACE with short jump FASTQ files trimmed with **[Sickle](https://github.com/najoshi/sickle) (v.1.33)** and reassess the goodness. In this way, we can determine whether to use the raw files (which have more data) or the trimmed files (which have better data) moving forward.  

Our last attempt at improving the assembly is closing the gaps between contigs using **GapCloser (v.1.12)** from SOAPdenovo2. We compare our de novo assembly with the known *Staphylococcus aureus* assembly using QUAST as a final step.  

## Results  

### Raw Data Analysis  
We started our analysis by looking at raw files and running FastQC. The short-insert files had 1557594 reads each, and the long-insert files had 1111432 reads each. Figures 1,2,3,4 indicate the per base quality assessment for raw files short insert size 1 and 2 and long insert size 1 and 2, respectively.  

#### Figure 1. FastQC Per Base Quality, Short-Insert 1  
![1](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/frag_1_fastqc.PNG)

#### Figure 2. FastQC Per Base Quality, Short-Insert 2
![2](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/frag_2_fastqc.PNG)  

#### Figure 3. FastQC Per Base Quality, Long-Insert 1  
![3](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/short_jump_1_fastqc.PNG)

#### Figure 4. FastQC Per Base Quality, Long-Insert 2  
![4](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/short_jump_2_fastqc.PNG)  

### Kmer Histogram Plots Before and After Correction  

We then used Jellyfish to evaluate sequencing error prevalence. We plot the number of times a 31-mer occurs on the x-axis, and number of 31-mers for each count on the y-axis for the short-insert library. The histogram plot before correction is shown in Figure 5.  

#### Figure 5. Jellyfish Histogram Before Correction  
![5](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/frag_1_31.pdf)  

The bin that corresponds to the initial valley point is bin 5, with 34039 kmers fitting that category.  

After correction of the short-insert files with SOAPdenovo2 modules, we replot the histogram, shown in Figure 6.  

#### Figure 6. Jellyfish Histogram After Correction  
![6](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/frag_1_31corrected.pdf)  

We go on to evaluate the corrected data files. Both short-insert corrected files had 1054393 reads. The new valley point is the bin 2, with frequency 15515. The peak of the graph is the bin 14, with frequency 142225. The corrected files 1 and 2 had average read length of 95.631 and 93.333 respectively, and total number of base pairs 100832704 and 98410110 respectively. Evaluating the genome size using the first file, we estimate it is 5244177 base pairs, approximately.  

**Extra Credit** The real genome size is around 3 million bp, so we are not very close to the the correct estimate. This is likely due to either an overestimate of the total bases or an underestimate of the depth of coverage. To determine the issue we would look into the correct depth of coverage and the mechanism of correction for kmer peak.  

### Assembly Goodness After Each Step  

In Table 7, we show how we evaluated the effectiveness of our assembly after each step of the process. We performed gap filling on the untrimmed short-read files, because the N50 score was higher before trimming.  

#### Table 7. Effectiveness of Each Step  

step|number of contigs|maximum contig length|N50  
---|---|---|---  
Contig assembly|699|80334|25284  
After Scaffolding Without Trimming|284|258872|128380  
After Scaffolding With Trimming|414|120417|41285  
After Gap Filling on Untrimmed|284|259114|128193  

Using class data, we plotted max contig length versus kmer length in Figure 8, and N50 score versus kmer length in Figure 9.  

#### Figure 8. Max contig length versus Kmer length  
![8](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/max_contig_vs_kmer_length.PNG)  

#### Figure 9. N50 versus Kmer length  
![9](https://github.com/cse185-sp18/cse185-week3-aarthivenkat/blob/master/labreport/n50_vs_kmer_length.PNG)  

### QUAST Evaluation  

Finally, we answer questions about the final QUAST report which compares our de novo assembly to the reference.  

**Why are there two results — what is the broken result?**  
According to the [QUAST manual](http://quast.bioinf.spbau.ru/manual.html#faq_q6) page, the "broken assembly" referes to QUAST splitting the fasta by contiguous fragments of Ns, and thus recontructing contigs used for construction of the scaffolds. Comparing the output and the output_broken can help the user determine if scaffolding was useful.  

**What fraction of the genome was covered by our scaffolds?**  
Comparing the total length of the output versus output_broken, 1.6% of the genome is covered by scaffolds.  

**What is the mismatch rate?**  
Considering the number of mismatches in the aligned bp versus the total aligned length, the mismatch rate is 0.0064%.  

**How many misassemblies are there?**  
There are 0 global misassemblies and 39 local misassemblies.  

**What are the orange/red scaffolds in the Icarus browser?**  
These refer to misassembled blocks.  


## Discussion  

In order to evaluate the success of this assembly, we delve further into the QUAST report. The genome fraction, or the ratio of aligned bases to the reference genome, is 97.688%. However, there were 39 local misassembly blocks, which, as evident through the contig browser, make up a large portion of our genome. At this high level, it is clear this is not a finished genome.  

QUAST predicted 2726 genes for my de novo genome using GeneMarkS, GeneMark-ES, MetaGeneMark, or GlimmerHMM, and 2846 using the inputted gene list. There are 2951 genes in the actual reference. Of those 2846, 49 are partially contained.  

To improve the assembly, we should begin with short-insert size fragment files with better per base quality scores, because our attempts to improve the quality with trimming failed. This could be achieved through higher coverage procedures. Further, we can work with different gap closing software tools to try to improve the N50 score through gap closing, because the N50 score went down after we used GapCloser on the raw short read files. Another approach would be to work with long read technologies, such as from Nanopore, in conjunction with short read technologies to leverage benefits in both quality and price. Finally, we could design libraries with varying insert sizes to improve scaffolding. 

Excluding some outliers due to errors in my classmates’ data, the general trend seemed to be a steady increase in both max contig length and N50 score as the kmer size increases. It seems that as we increase kmer size, we have a more acute error correction procedure, but we are still limited by technology, as increasing kmer size also makes de Brujin graphs more difficult to create.  

## Citations  

1. Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi and Glenn Tesler, QUAST: quality assessment tool for genome assemblies, Bioinformatics (2013) 29 (8): 1072-1075.doi: 10.1093/bioinformatics/btt086  
2. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc  
3. Guillaume Marcais and Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics (2011) 27(6): 764-770  
4. Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files (Version 1.33) [Software]. Available at https://github.com/najoshi/sickle.  
5. Luo, Ruibang et al. “SOAPdenovo2: An Empirically Improved Memory-Efficient Short-Read de Novo Assembler.” GigaScience 1 (2012): 18. PMC. Web. 23 Apr. 2018.  
6. K. Salikhov, G. Sacomoto and G. Kucherov. Using cascading Bloom filters to improve the memory usage for de Brujin graphs, WABI 2013  
7. Marten Boetzer, Christiaan V. Henkel, Hans J. Jansen, Derek Butler, Walter Pirovano; Scaffolding pre-assembled contigs using SSPACE, Bioinformatics, Volume 27, Issue 4, 15 February 2011, Pages 578–579, https://doi.org/10.1093/bioinformatics/btq683  
8. Nagarajan, N. and Pop, M (2013). "Sequence assembly demystified" Nature reviews. Genetics 14(3): 157-167.  
