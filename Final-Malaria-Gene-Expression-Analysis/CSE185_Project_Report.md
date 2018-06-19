# Whole Transcriptome Differential Expression Analysis of Chloroquine-Resistance Mechanisms in *Plasmodium Falciparum*
#### Aarthi Venkat | 10 June 2018

## Abstract
Despite decades of research of malaria, there is still much left unknown about malarial parasites and their modes of resistance. In this paper, I attempt to employ NGS strategies to add to the literature of the mechanism of chloroquine resistance. Through alignment-independent measures, I quantified RNA-Seq reads of two resistant strains of *Plasmodium falciparum*, one cultured in the presence of chloroquine and one in the absence. I then identified 4784 differentially expressed transcripts and characterized their ontology as highly associated with movement/invasion as well as glucose metabolism. Through my analysis, I not only reproduced the results of the paper from which I obtained the reads with a different pipeline, but also extended the discussion and analysis of gene ontology.

## Introduction  
According to the World Health Organization[<sup>[**1**]</sup>](http://www.who.int/malaria/media/world-malaria-report-2017/en/), 216 million cases of malaria occurred worldwide in 2016, with 445,000 of these resulting in death. Among all malarial parasites, *Plasmodium falciparum* is the most insidious, accounting for 99% of estimated malaria cases in 2016. Therefore, it is a critical public health concern to understand the mechanism of the malarial parasite as well as subsequent mechanisms of resistance, particularly of *P. falciparum*, which is known to develop resistance against common treatments.  

Researchers have learned that inside red blood cells, the malarial parasite must degrade hemoglobin to acquire amino acids for its own protein construction and metabolism. During this process, the parasite releases a toxic molecule heme, and, to avoid self-destruction, crystallizes it to form hemozoin, a nontoxic molecule[<sup>[**2**]</sup>](https://www.ncbi.nlm.nih.gov/pubmed/11915943).  

Hemozoin, thus, is an integral part of the malarial mechanism and is a desirable target to treat malaria[<sup>[**3**]</sup>](https://www.sciencedirect.com/science/article/pii/S0304416514000555?via%3Dihub). Enter chloroquine, a prevalent anti-malarial medication, which diffuses into the red blood cell and caps hemozoin molecules to prevent further biocrystallization of heme, leading to heme buildup. Chloroquine then binds to heme to form a highly toxic complex which disrupts membrane function and ultimately results in cell lysis and autodigestion[<sup>[**4**]</sup>](https://www.sciencedirect.com/science/article/pii/016372589390056J?via%3Dihub). Chloroquine has been extensively used to treat malaria, though *P. falciparum* has begun to develop resistance to the treatment<sup>[**4**]</sup>. Thus, it is of high interest to analyze *P. falciparum*'s mode of resistance, and recent emphasis has been placed on NGS strategies to better characterize the mechanism.  

A recent study from Antony et al performed a whole transcriptome expression analysis with RNA-Seq between three samples of *P. falciparum*, one chloroquine sensitive (3D7), one chloroquine resistant exposed to chloroquine (Dd2_CQ), and one chloroquine resistant not exposed to chloroquine (Dd2)[<sup>[**5**]</sup>](https://www.sciencedirect.com/science/article/pii/S2213596016300496). In this paper, I will use the RNA-Seq data acquired from Antony et al and perform transcript quantification and differential expression analysis using an alignment-independent pipeline (Kallisto-Sleuth), contrasting to the alignment-dependent pipeline (TopHat-Cufflinks-Cuffdiff) of the Antony et al paper. At the end of my analysis, I will compare my results of Dd2 vs. Dd2_CQ expression with that of the Cuffdiff output from the paper, as well as extend the analysis to further characterize gene ontology from differentially expressed genomic regions. Through this, I hope to add to the literature of the molecular mechanism of malaria drug resistance and contribute towards better understanding of potential drug targets.  

## Methods

### RNA-Seq Datasets
The data from the original paper is paired-end RNA-Seq data from Illumina HiSeq 2500 sequencing. Reads were collected from three samples of Plasmodium falciparum - one non-resistant strain (Pfal3D7;SRX1546678) and two resistant, one cultured in the presence of chloroquine (PfalDd2_CQ;SRX1535489) and one in the absence (PfalDd2;SRX1546677). I also procured the differential gene expression results from GEO (GSE77499) in the form of the output of Cuffdiff, so that I could perform downstream gene ontology analyses to compare with my pipeline's output. 

For the sake of this project, I worked only with the resistant strain, comparing reads from PfalDd2 and PfalDd2_CQ (referred to as Dd2 and Dd2_CQ in this paper). I downsampled the data to the first 1.5 million reads of each paired-end file. I then ran [**FastQC (v0.11.7)**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<sup>[**6**]</sup> on every FASTQ file. 

### Transcript Quantification  
The next step of my analysis involved transcription quantification through an alignment-independent pipeline. I opted to use [**Kallisto (v0.44.0)**](https://pachterlab.github.io/kallisto/about)<sup>[**7**]</sup>, first indexing with default parameters the *Plasmodium falciparum* transcriptome obtained from [**PlasmoDB (Release 37)**](http://plasmodb.org/plasmo/)<sup>[**8**]</sup> through this [script](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/scripts/fasta.sh). Then I ran kallisto quant with the index, 3 threads and 100 bootstrap samples. I invoked multithreading to increase the speed, and used 100 bootstrap samples as a way to sample small sets from a larger dataset and develop reliable statistics about transcript abundance.  

As an additional test, I compared the percentage of reads pseudoaligned and compared it to the percentage of reads aligned with Tophat in the original paper's analysis. I then performed pairwise correlation analysis using a custom [script](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/scripts/correlation.sh). In the custom script, I pasted together the abudance.tsv output from kallisto, extracted just the tpm values, removed the "tpm" header, extracted only values where both columns are non-zero, and ran a pearson correlation test.  

### Differential Expression Analysis  
Following Kallisto, I used a [script](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/scripts/sleuth) invoking [**Sleuth**](https://pachterlab.github.io/sleuth/about)<sup>[**9**]</sup>, which interprets the kallisto output, writes the significant differentially expressed transcripts to a tab-delimited file. The program required I use replicates, so I created pseudoreplicates by copying the sample as a separate library. The False Discovery Threshold is 0.05. I then determined the number of significant transcripts and created several visualizations of transcript abundance through the sleuth_live interactive feature.  

### Gene Ontology  
To better functionally interpret the significant transcripts, I turned to the [**PlasmoDB** Gene Identification Feature](http://plasmodb.org/plasmo/showQuestion.do?questionFullName=GeneQuestions.GeneByLocusTag), where I searched the top 10 transcript IDs from my output, resulting in a summary list of the genes and gene ontology terms. I then invoked [**REVIGO (Jan 2017 release)**](http://revigo.irb.hr/)<sup>[**10**]</sup>, a new tool to the class which reduces the gene ontology list by summarizing redundant terms and visualizing the list in a word cloud to better understand function.  

### Comparison to Original Paper through Gene Ontology 
In order to compare my results to the downstream results from the original paper, I ran a [script](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/scripts/downsample.sh) to downsample the original paper's Cuffdiff output to just Dd2 and Dd2_CQ comparisons with significance at the same threshold, q<=0.05. I then determined the number of transcripts, and put all the transcripts into PlasmoDB to study gene ontology. When that query returned few results, I put the list of genes the Cuffdiff output implicated into PlasmoDB, and returned a summary list of the genes and gene ontology terms as before. Finally, I used REVIGO on the GO terms to make a word cloud, so that I could compare the word clouds between my analysis and the paper's.

## Results  

### Quality Check  
After running FastQC, I was able to produce the following plots in **Figure 1** for sequency quality and sequence content for each file in each paired-end library. It seems that the reads are all with quality >= 30 Phred score, and the adapters are cut, which leads me to believe the data is already processed through quality control. Furthermore, every file failed only one module: **Per Base Sequence Content**. This is because the difference between A and T or G and C was greater than 20% within the first 12 bases. This error was elicited by biased fragmentation, due to the generation of the libraries from ligation of random hexamers. These libraries thus have an intrinsic selection bias in the positions at which reads start. According to the [**FastQC documentation**](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html), this does not represent biased individual sequences, and while nearly all RNA-Seq libraries will fail this module, it is not a problem that can be fixed through processsing, and it does not affect downstream analysis.  

#### Figure 1. FastQC Output  
**Dd2_1**  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2_1_quality.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2_1_content.PNG)  
**Dd2_2**  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2_2_quality.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2_1_content.PNG)  
**Dd2_CQ_1**  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2CQ_1_quality.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2CQ_1_content.PNG)  
**Dd2_CQ_2**  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2CQ_1_quality.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/fastqc/Dd2CQ_1_content.PNG)  

### Kallisto Alignment & Correlation Tests  
After running Kallisto Quant to quantify transcript abundance, I could calculate the number of reads pseudoaligned and compare that to the original paper's alignment statistics in **Table 2**. The alignments rates are comparable and high.

#### Table 2. Percentage of Reads Aligned  
Sample|My Results|Paper Results  
---|---|---  
Dd2|91.28%|89.19%  
Dd2_CQ|89.19%|84%  

The concordance of the two samples is **0.483182958**.  

### Sleuth Results  
After running Sleuth, the number of significantly differentially expressed transcripts was **4784**, whereas for the original paper it was **29** at the same significance threshold. Through the sleuth_live interactive visualization feature, I could visualize a table of all the differentially expressed transcripts, shown in **Table 3**. After identifying highly differential transcripts, I visualized box plots of abundance differences for the first three in **Figure 4**, as well as a heatmap of the first 10 transcripts in **Figure 5**. Finally, I took advantage of the interface to look at summaries of Kallisto output, namely the abundance distributions and fragment length distributions for the two samples (**Figure 6**, **Figure 7**, respectively).  

#### Table 3. Table of Differentially Expressed Transcripts  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/table.PNG)  

#### Figure 4. Box Plots of Abundance for Top Three DE Transcripts  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/transcript_abundances_1.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/transcript_abundances_2.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/transcript_abundances_3.PNG)  

#### Figure 5. Heatmap of Top Ten DE Transcripts  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/transcript_abundances_heatmap.PNG)  

#### Figure 6. Abundance Distributions from Two Samples  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/abundance_distribution.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/frag_length_distribution_87.PNG)  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/sleuth/frag_length_distribution_88.PNG)  

### Gene Ontology  
From PlasmoDB, I was able to get a summarized list of the top ten transcripts and their descriptions (**Table 7**), as well as perform a GO analysis and get a summarized list of the gene ontology (**Table 8**). I could then utilize REVIGO to form a word cloud from the GO IDs (**Figure 9**).  

#### Table 7. Transcript Descriptions  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/my_gene_list.PNG)  

#### Table 8. Gene Ontology  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/my_gene_ontology.PNG)  

#### Figure 9. Word Cloud  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/my_revigo_wordcloud.PNG)  

In order to compare to original paper's results, I performed a similar analysis with their differentially expressed transcripts output, located from the GEO. I first used their 29 significantly DE transcripts as input, getting three hits from PlasmoDB (**Table 10**). I then used the genes associated with those transcritps as input, getting seven hits this time (**Table 11**). From this second list, I performed GO analysis and got a summarized list of the gene ontology (**Table 12**). Finally, I again used REVIGO to form a word cloud from the GO IDs (**Figure 13**).  

#### Table 10. Original Paper Transcript Descriptions - Transcript ID Input  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/original_gene_list.PNG)  

#### Table 11. Original Paper Transcript Descriptions - Gene ID Input  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/original_gene_list_geneidsearch.PNG)  

#### Table 12. Original Paper Gene Ontology  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/original_gene_ontology.PNG)  

#### Figure 13. Original Paper Word Cloud  
![](https://github.com/cse185-sp18/cse185-final-project-aarthivenkat/blob/master/results/GO/original_revigo_wordcloud.PNG)  

## Discussion 
### Analysis of Gene Ontology from Kallisto-Sleuth Output  
The significant gene ontology results suggest **entry/movement into the host cell** to be the top function of differentially expressed genes, among immune system and metabolic processes. Similarly, the word cloud identifies **invasion**, **entry**, **penetration**, and **growth** among keywords correlating with the GO terms. 

Knowing that chloroquine diffuses into the membrane of the digestive vacuole of the parasitic cell and forms a toxic complex which disrupts membrane function, it is clear that the terms involving invading or moving are associated with the differential membrane function; this makes sense and is a good sanity check on our results.    

Another metabolic process that seems to stand out is glucose metabolism, namely from the GO terms **gluconeogenesis**, **hexose biosynthetic/metabolic process**, **carbohydrate biosynthetic process**, and **monosaccharide biosynthetic/metabolic process**. This connection was not quite as intuitive to me, but upon research it seems that the effect of chloroquine on insulin and glucose homeostasis is a well-documented phenomenon, dating back to the 1980s, where researchers found chloroquine-intake associated with decreased degradation of insulin[<sup>[**11**]</sup>](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1245517/), and significant changes in glucose tolerance[<sup>[**12**]</sup>](https://www.ncbi.nlm.nih.gov/pubmed/12087927). Thus, it seems that when the parasite is exposed to chloroquine, even if, through the resistance mechanisms, the cell doesn't lyse, the chloroquine still affects the carbohydrate metabolic processes.  

### Comparison to Original Paper
Both my paper and the original paper had a high percentage of reads aligned. However, my analysis outputted 4784 transcripts, while the original anlsyis outputted 29. This discrepency could be due to the difference in the pipeline or the downsampling of reads in my dataset. Of the 29 transcripts IDs inputted into PlasmoDB from the Cuffdiff output, 3 were recognized by the database, and 2 were in my list of the top ten differentially expressed transcripts form Sleuth. Thus, there exists some indication of overlap with results between the alignment-dependent TopHat-Cufflinks-CuffDiff approach, and the alignment-independent Kallisto-Sleuth approach. 

As for comparison for gene ontology - my analysis included largely both ontology terms related to invasion as well as glucose metabolism. However, the original paper had ontology terms only related to invasion. The word clouds further exemplify this difference, with my word cloud containing a whole host of terms that could relate to function beyond invasion, such as **nucleoside** and **monosaccharides**, while the original paper's word cloud only including versions of **parasitism** and **locomotion**. The reason for this difference is potentially that with many more transcripts determined as "differentially expressed", I had a wider variety of transcripts with different functions and "locomotion" did not dominate the sphere of gene ontology as it did for the original paper output. Ultimately, I find it interesting that my analysis was able to still identify key gene functions associated with chloroquine, beyond even the original paper's output, proving that my RNA-Seq analysis was accurate in many respects and could be implemented as a pipeline for future projects.  

### Further Research
It is necessary to discuss that the concordance between two samples that are both PfalDd2 was extremely low (0.483), when I was expecting a relatively high concordance (>0.8). Further, I received nearly 5000 differentially expressed transcripts, which is a very large number considering the anticipated similarity between the two genomes. These discrepancies suggest that there was a lot of noise in my data, and with only 1.5 million reads and no replicates, Kallisto and Sleuth were not able to sufficiently filter out the noise. Thus, for further research, RNA-Seq data requires technical replicates and biological replicates in order to improve the accuracy of the tools being used. It would also be interesting to integrate time-dependent data to see how the genes are expressed as the organism is exposed to chloroquine for an extended period of time. Finally, adding analyses such as protein conformation, chromatin conformation, and other epigenetic studies could truly extend the reach of *Plasmodium falciparum* and chloroquine-resistance research.

## Citations
1. “Key Points: World Malaria Report 2017.” World Health Organization, World Health Organization, 14 Dec. 2017, www.who.int/malaria/media/world-malaria-report-2017/en/.  
2. Beeson, J G, and G V Brown. “Pathogenesis of Plasmodium Falciparum Malaria: the Roles of Parasite Adhesion and Antigenic Variation.” Advances in Pediatrics., U.S. National Library of Medicine, Feb. 2002.  
3. Coronado, Lorena M., Christopher T. Nadovich, and Carmenza Spadafora. “Malarial Hemozoin: From Target to Tool.” Biochimica et biophysica acta 1840.6 (2014): 2032–2041. PMC. Web. 11 June 2018.  
4. Slater, Andrew. “Chloroquine: Mechanism of Drug Action and Resistance in Plasmodium Falciparum.” Egyptian Journal of Medical Human Genetics, Elsevier, 14 Nov. 2002.  
5. Antony, Hiasindh. “Whole Transcriptome Expression Analysis and Comparison of Two Different Strains of Plasmodium Falciparum Using RNA-Seq.” Egyptian Journal of Medical Human Genetics, Elsevier, 23 Apr. 2016.
6. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc  
7. Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519  
8. PlasmoDB: a functional genomic database for malaria parasites. Nucleic Acids Res. 2008 Oct 31. Aurrecoechea C, et al.  
9. Harold J. Pimentel, Nicolas Bray, Suzette Puente, Páll Melsted and Lior Pachter, Differential analysis of RNA-Seq incorporating quantification uncertainty, Nature Methods (2017), advanced access http://dx.doi.org/10.1038/nmeth.4324.  
10. Supek F, Bošnjak M, Škunca N, Šmuc T. "REVIGO summarizes and visualizes long lists of Gene Ontology terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800  
11. Smith, G D et al. “Effect of Chloroquine on Insulin and Glucose Homoeostasis in Normal Subjects and Patients with Non-Insulin-Dependent Diabetes Mellitus.” British Medical Journal (Clinical research ed.) 294.6570 (1987): 465–467. Print.  
12. Gaafar, Khadiga M. “Effect of Chloroquine on Glucose Metabolism.” Arzneimittelforschung, Thieme Medical Publishers, 26 Dec. 2011.  
