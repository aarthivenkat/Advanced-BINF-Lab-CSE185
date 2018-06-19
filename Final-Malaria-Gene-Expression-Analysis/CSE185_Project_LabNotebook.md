# Lab notebook for CSE185 Final Project

## Week 8 Tuesday
**3 relevant datasets**

I'm interested in researching malaria through differential gene expression for mosquito strains and Plasmodium parasitic strains. I am currently looking at papers involved in functional analysis, PCA, and other RNA-Seq analysis components. 

1. Differential gene expression and PCA - RNA-Seq pipelining of mosquito strain  
https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-018-2817-5#Abs1  
SRA RNA-Seq reads: PRJNA390156  
RNA-Seq pipeline: https://github.com/bjmain/Acol_cyp1_DESeq2/blob/master/DE_analysis_supp_file_bmain.ipynb  

2. Whole Transcriptome Expression Analysis of Plasmodium Falciparum  
https://www.sciencedirect.com/science/article/pii/S2213596016300496  
Raw data: BioProject: 308455  
Processed data: GEO: GSE77499  

3. Functional Analysis of P.vivax Transcriptome  
https://www.nature.com/articles/srep20498  
GEO: GSE61252  


## Week 8 Thursday

I decided to choose dataset 2 because I had access to the raw data as well as the processed data, and the paper was very clear about the RNA-Seq pipeline, and I felt as though there was room for improvement in terms of the statistical analysis, as the paper only created volcano plots and did not attempt to characterize differential gene expression any further.   

The paper's method of processing the FASTQ reads involved:  
1. QC check with FastQC  
2. Base trimming through in-house Perl script (I could use Sickle or Trimmomatic)  
3. Adapter trimming with Cutadapt  
4. Contamination removal with bowtie2, various in-house scripts, and picardtools  
5. Read alignment with Tophat and default parameters  
6. Expression estimation with cufflinks  
7. Differential expression analysis with cuffdiff.  

As the paper did not go into detail regarding the contamination removal process, I find it necessary to not include that step in my pipeline, which would essentially be a replication of this. I will extend the pipeline and the analysis but incorporating **cummeRbund**, which I have installed locally and works with cuffdiff output to visualize and analyze comparative analyses.  

Important Links:  
Paper: file:///C:/Users/Aarthi/Downloads/1-s2.0-S2213596016300496-main.pdf  
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77499  
Pfal3D7 SRR: https://www.ncbi.nlm.nih.gov/sra/SRX1546678[accn]  
PfalDd2 SRR: https://www.ncbi.nlm.nih.gov/sra/SRX1546677[accn]  
PfalDd2_CQ SRR: https://www.ncbi.nlm.nih.gov/sra/SRX1535489[accn]  
Pipeline description of RNA-Seq: https://darwinawardwinner.github.io/resume/examples/Salomon/Teaching/RNA-Seq%20Lecture.pdf  
Pipeline analysis of RNA-Seq: https://darwinawardwinner.github.io/resume/examples/Salomon/Teaching/RNA-Seq%20Lab.html  
cummeRbund documentation: https://www.bioconductor.org/packages/devel/bioc/vignettes/cummeRbund/inst/doc/cummeRbund-example-workflow.pdf  

## Week 9 Tuesday

In order to resrict the memory usage, I've decided to only work with two datasets - PfalDd2 and PfalDd2_CQ (resistant strain and resistant strain exposed to chloroquine). I've taken the first 1.5 million reads of each one.  

```
fastq-dump --split-files -X 1500000 SRR310787 #Dd2_CQ  
fastq-dump --split-files -X 1500000 SRR310788 #Dd2  
```

Next, I run FastQC on the four files (two paired-end per dataset).  

```
fastqc -o ./FastQC-Output SRR3107187_1.fastq SRR3107187_2.fastq SRR3107188_1.fastq SRR3107188_2.fastq  
scp avenkat@ieng6.ucsd.edu:/home/linux/ieng6/cs185s/avenkat/final/*.html .  
```

Every file failed only one module: Per Base Sequence Content. This is because the difference between A and T or G and C was greater than 20% within the first 12 bases. This error was elicited by biased fragmentation, due to the generation of the libraries from ligation of random hexamers. These libraries thus have an intrinsic selection bias in the positions at which reads start. According to the ![FastQC documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html), this does not represent biased individual sequences, and while nearly all RNA-Seq libraries will fail this module, it is not a problem that can be fixed through processsing, and it does not affect downstream analysis.

## Week 9 Thursday

Today was spent installing TopHat and strengthening my pipeline using industry standards exemplified [here](https://www.illumina.com/documents/products/technotes/RNASeqAnalysisTopHat.pdf)

## Week 10 Tuesday

Realizing that I need an alignment-independent pipeline for time and memory purposes, I'm turning to kallisto.  

### NEW PIPELINE
PlasmoDB + FastQC + Kallisto + Sleuth

```
kallisto index -i transcripts.idx PlasmoDB-37_Pfalciparum3D7_AnnotatedTranscripts.fasta 
```
[build] loading fasta file PlasmoDB-37_Pfalciparum3D7_AnnotatedTranscripts.fasta   
[build] k-mer length: 31   
[build] warning: clipped off poly-A tail (longer than 10)   
        from 1 target sequences   
[build] counting k-mers ... done.  
[build] building target de Bruijn graph ...  done   
[build] creating equivalence classes ...  done  
[build] target de Bruijn graph has 52070 contigs and contains 11940447 k-mers   

```
for i in SRR3107187 SRR3107188; do kallisto quant -t 3 -b 100 -o $i -i transcripts.idx ${i}_1.fastq ${i}_2.fastq; done
```
[quant] fragment length distribution will be estimated from the data   
[index] k-mer length: 31   
[index] number of targets: 5,800   
[index] number of k-mers: 11,940,447  
[index] number of equivalence classes: 18,600  
[quant] running in paired-end mode  
[quant] will process pair 1: SRR3107187_1.fastq  
                             SRR3107187_2.fastq  
[quant] finding pseudoalignments for the reads ... done  
[quant] processed 1,500,000 reads, 1,369,135 reads pseudoaligned  
[quant] estimated average fragment length: 165.954  
[   em] quantifying the abundances ... done  
[   em] the Expectation-Maximization algorithm ran for 851 rounds  
[bstrp] number of EM bootstraps complete: 100  


[quant] fragment length distribution will be estimated from the data  
[index] k-mer length: 31  
[index] number of targets: 5,800  
[index] number of k-mers: 11,940,447  
[index] number of equivalence classes: 18,600  
[quant] running in paired-end mode  
[quant] will process pair 1: SRR3107188_1.fastq  
                             SRR3107188_2.fastq  
[quant] finding pseudoalignments for the reads ... done  
[quant] processed 1,500,000 reads, 1,303,785 reads pseudoaligned  
[quant] estimated average fragment length: 166.382  
[   em] quantifying the abundances ... done  
[   em] the Expectation-Maximization algorithm ran for 547 rounds  
[bstrp] number of EM bootstraps complete: 100  

```
paste SRR3107187/abundance.tsv SRR3107188/abundance.tsv | cut -f 5,10 | grep -v tpm | awk '(!($1==0 && $2==0))' | datamash ppearson 1:2
```
0.4831829580023   
It seems that the samples are not very concordant!  

## Week 10 Thursday  

Sleuth took years to get working, I had to create pseudoreplicates in order for it to run - today I learned many RNA-Seq tools require replicates in order to garner accurate results!

```
> library("sleuth")
> sample_id = c("SRR3107187", "SRR3107187_1", "SRR3107188", "SRR3107188_1")
> kal_dirs = file.path(sample_id)
> s2c = read.table(file.path("C:/Users/Aarthi/Documents/exp_info.txt"), header = TRUE, stringsAsFactors=FALSE)
> s2c = dplyr::mutate(s2c, path = kal_dirs)
> so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
reading in kallisto results
dropping unused factor levels
....
normalizing est_counts
5072 targets passed the filter
normalizing tpm
merging in metadata
summarizing bootstraps
....
Warning message:
In check_num_cores(num_cores) :
  It appears that you are running Sleuth from within Rstudio.
Because of concerns with forking processes from a GUI, 'num_cores' is being set to 1.
If you wish to take advantage of multiple cores, please consider running sleuth from the command line.
> so = sleuth_fit(so, ~condition, 'full')
fitting measurement error models
shrinkage estimation
computing variance of betas
> so = sleuth_fit(so, ~1, 'reduced')
fitting measurement error models
shrinkage estimation
2 NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit.
The LOESS fit will be repeated using exact computation of the fitted surface to extrapolate the missing values.
These are the target ids with NA values: PF3D7_0532000.1, PF3D7_0726000.1
computing variance of betas
> so = sleuth_lrt(so, 'reduced', 'full')
> sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
> sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
> write.table(sleuth_significant, "C:/Users/Aarthi/Documents/sleuth_results.tab",  sep="\t", quote=FALSE)
```
To determine how many significant transcripts there are,
`wc -l sleuth_results.tab` - 1 for header = 4784  

While this is a lot of transcripts, I don't know if that's a good thing. The fact that there are no replicates means any small difference is likely significant. Nonetheless, I want to use the capabailities of sleuth_live and RShiny to explore this data further.  
I've done this, taking screenshots of the results & plotting the top transcripts.

I want to see Volcano plots, so I run the Wald test.
```
so <- sleuth_wt(so, 'conditionDd2_CQ')
sleuth_live(so)
```
Next stop to figure out the genes and ontology associated with my results - PlasmoDB!  

## Outside Lab  - PlasmoDB for Gene Ontology

Today I learned - PlasmoDB is a godsend. I was able to get a gene list and gene ontology results with ease. I even could create a word cloud!

## Outside Lab - Start analyzing paper results

Now, I have collated a whole host of information about my results, but I have to find statistically sound ways to compare the results. 

(1) Since I downsampled the results to only Dd2 and Dd2_CQ, I have to downsample the differential gene expression results as well. I set the significance threshold as the same as for sleuth = q-val <= 0.05.  
 ```
 cat GSE77499_Dd2CQ_Dd2_3D7_gene_exp.diff | awk '($5=="Dd2CQ" && $6=="Dd2" && $13<=0.05)' > original_only_dd2.diff
 ```
 `wc -l original_only_dd2.diff` 29 differentially expressed transcripts.
 
 I did gene ontology for the original, and the two rRNA transcripts were also in my analysis. Huzzah! There was also an snRNA transcript I did not catch.  
 
 The gene ontology tool did not work for these transcripts ("Your result has no genes with GO terms, so you can't use this tool on this result. Please revise your search and try again."). Bummer. PlasmoDB for the genes from the Cuffdiff file?  
 
 I need one more tool. I could throw some IGV in there, or I could use something similar that is definitely a new tool.
 REVIGO Word Cloud!
 
## Data Check

Data in ieng6 at the end of it all (**1.8 GB**)
drwxr-s--- 2 avenkat ieng6_cs185s 4.0K Jun  9 15:33 --gtf
-rw-r----- 1 avenkat ieng6_cs185s 2.0M Jun  9 20:17 GSE77499_Dd2CQ_Dd2_3D7_isoform_exp.diff
-rw-r----- 1 avenkat ieng6_cs185s  15M Apr 19 12:36 PlasmoDB-37_Pfalciparum3D7_AnnotatedTranscripts.fasta
-rw-r----- 1 avenkat ieng6_cs185s 397M May 31 12:31 SRR3107187_1.fastq
-rw-r----- 1 avenkat ieng6_cs185s 397M May 31 12:31 SRR3107187_2.fastq
-rw-r----- 1 avenkat ieng6_cs185s 397M May 31 12:30 SRR3107188_1.fastq
-rw-r----- 1 avenkat ieng6_cs185s 397M May 31 12:30 SRR3107188_2.fastq
-rw-r----- 1 avenkat ieng6_cs185s 3.4K Jun  9 20:33 original_only_dd2.diff
-rw-r----- 1 avenkat ieng6_cs185s 244M Jun  9 15:27 transcripts.idx
