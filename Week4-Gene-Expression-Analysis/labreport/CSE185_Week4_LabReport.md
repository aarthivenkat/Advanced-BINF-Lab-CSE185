# Deletions in Snake ZRS enhancer to *Shh* gene implicated in lack of snake limb development  
##### Aarthi Venkat
##### Date: May 1 2018

## Abstract  

We studied hindlimb, forelimb, and midbrain RNA-Seq mice samples to determine which genomic regions could be implicated in limb development. Through transcript quantification and analysis, we determined 1940 differentially-expressed genes and focused on *Shh*, critical in early development. To analyze regulatory regions for *Shh*, we extend our analysis to ChIP-Seq and PhyloP data, identifying the ZRS region as highly conserved and involved in limb-specific regulation. Finally, we compare sequence alignment of ZRS across species, implicating two deletion regions in exclusively snake sequences as factors in limb development. These findings are important in understanding embryonic development regulation and fin-to-limb evolutionary genomics.  

## Introduction  

In this lab, we are interested in determining which genes and which regulatory regions are implicated in limb development, as well as whether mutations in these critical regions could potentially explain the absence of limbs in pythons, rattlesnakes, cobras, and boas. To answer this question, we first observe hindlimb, forelimb, and midbrain RNA-Seq data from mouse early development by quantifying expression and comparing expression patterns in significant genes across datasets, and then we analyze ChIP-Seq and PhyloP data to integrate histone modification and conservation information, respectively. Ultimately, we reason that by performing differential analysis in limb tissues versus non-limb tissues, we can determine which regions may critical for limb development, and we can identify such regions and their expression in other limbed organisms as well as non-limbed organisms (i.e. snakes) to observe evolutionary differences.  

The question of which genomic regions are implicated in limb development is a critical one in terms of evolutionary biology. According to [Lettice et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5561167/)<sup>[**4**]</sup>, embryonic development is a highly spatiotemporally regulated process. To study early development is to study the structural role of enhancer activity, and studying evolution and disease requires studying early development. Hence, in order to make contentions about evolutionary biology, we must study enhancer architecture (such as the ZRS region for limb-specific expression) and critical development genes (such as the *Shh* gene). The conservation of the network of expression contributes to our understanding of the evolutionary transition from fin-to-limb and can answer comparative interspecies questions at large. Further, we are interested in this question so we can understand the genomic basis behind diseases arising in early development.  

## Methods  

### RNA-seq datasets  

I am working with chr5 RNA-seq data from 3 tissues and two replicates per tissue. The three tissues are from a developing mouse: hindlimb (HL), forelimb (FL), and midbrain (MB). Each replicate has two FASTQ files associated with it, indicating that the datasets are paired. I determined the number of reads and read length by running a [custom script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/custom_scripts) (line 1); the results are in **Table 1**. After inspecting the files, I ran [**FastQC (v0.11.7)**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<sup>[**1**]</sup> on every fq.gz file using another [custom script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/custom_scripts) (line 2). Every file failed only one module: **Per Base Sequence Content**. This is because the difference between A and T or G and C was greater than 20% within the first 12 bases. This error was elicited by biased fragmentation, due to the generation of the libraries from ligation of random hexamers. These libraries thus have an intrinsic selection bias in the positions at which reads start. According to the [**FastQC documentation**](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html), this does not represent biased individual sequences, and while nearly all RNA-Seq libraries will fail this module, it is not a problem that can be fixed through processsing, and it does not affect downstream analysis. Finally, from the BAMfile header of reads aligned for us, we can interpret using a [custom script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/custom_scripts) (line 3) that the reference genome build used throughout the lab is mm10 for Mouse.  

#### Table 1. Number of Reads and Read Length for each Dataset  
File|NumReads|ReadLength  
---|---|---  
FL_Rep1_chr5_1.fq.gz|3179399|50  
FL_Rep1_chr5_2.fq.gz|3179399|50  
FL_Rep2_chr5_1.fq.gz|2965029|50  
FL_Rep2_chr5_2.fq.gz|2965029|50  
HL_Rep1_chr5_1.fq.gz|3932838|50  
HL_Rep1_chr5_2.fq.gz|3932838|50  
HL_Rep2_chr5_1.fq.gz|2811913|50  
HL_Rep2_chr5_2.fq.gz|2811913|50  
MB_Rep1_chr5_1.fq.gz|3254975|50  
MB_Rep1_chr5_2.fq.gz|3254975|50  
MB_Rep2_chr5_1.fq.gz|3413939|50  
MB_Rep2_chr5_2.fq.gz|3413939|50  

### RNA-seq analysis  

Now, we move towards transcript abundance quantification with [**Kallisto (v0.44.0)**](https://pachterlab.github.io/kallisto/about)<sup>[**6**]</sup>. Using a [custom script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/run_kallisto.sh), we run kallisto quant with a GTF file with transcriptome information and a kallisto index, as well as 3 threads and 100 bootstrap samples. We invoke multithreading to increase the speed, and we use 100 bootstrap samples as a way to sample small sets from a larger dataset and develop reliable statistics about transcript abundance. Then, we complete pairwise correlation analysis using a [custom script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/custom_scripts) (line 4). In the custom script, we paste together the abudance.tsv output from kallisto, extract just the tpm values, remove the "tpm" header, extract only values where both columns are non-zero, and run a pearson correlation test. Finally, we used a [script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/sleuth_script) invoking [**Sleuth**](https://pachterlab.github.io/sleuth/about)<sup>[**11**]</sup>, which interprets the kallisto output and writes the significant differentially expressed transcripts to a tab-delimited file. The False Discovery Threshold is 0.05.  

### Enhancer analyses    

On [**IGV (v2.4)**](http://software.broadinstitute.org/software/igv/)<sup>[**2**]</sup>, we loaded 6 given tdf tracks (2 replicates for each of 3 tissues), where the tdf files demonstrate the read counts per position for the entire genome. Then, we loaded 6 given bedGraph files from ChIP-seq experiments for H3K27ac ad H3K4me1, histone modifications found near enhancer regions. These bedGraph files will be useful in identifying putative enhancers in the regulatory regions of *Shh*. Finally, we add one additional track to IGV - the PhyloP track, which gives a per-base pair measure of sequence conservation. We got the PhyloP track from [**UCSC Genome Table Browser**](http://genome.ucsc.edu/cgi-bin/hgTables)<sup>[**9**]</sup>. Finally, we intend to perform multiple sequence alignment of the ZRS region across species, as we have deduced that this is an important region for limb development. To do so, we invoke [**MAFFT (v7.397)**](https://mafft.cbrc.jp/alignment/software/)<sup>[**3**]</sup> and [**MView (v1.63)**](https://www.ebi.ac.uk/Tools/msa/mview/)<sup>[**7**]</sup>. There were no non-default parameters for MAFFT, and the non-default parameters for MView were `-html full` to produce an HTML of the alignment, and `-coloring any` to color the nucleotides in the HTML.  

## Results  

### RNA-seq analysis  

The replicates were highly concordant, with correlation rate of 0.995. All replicates were more concordant with each other than with other tissues. Regarding comparisons between tissues, FL and HL had a correlation rate of 0.959, whereas FL and MB 0.954 and HL and MB 0.931. Therefore, FL and HL were the most similar, and HL and MB the most different. Figure 2 highlights the correlation analysis for each of the tissues and replicates.  

#### Figure 2. Pairwise Analysis Across Tissues and Replicates  
![correlation analysis](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/correlation.PNG)  

Regarding differentially expressed transcripts, I found 1940 significant transcripts. The top 10 genes and their names (collated from a given [**Ensembl mouse index**]( http://uswest.ensembl.org/Mus_musculus/Info/Index)<sup>[**8**]</sup>) are listed in Table 3 below. For the first top 5 hits, I navigated to that gene in IGV and took screenshots, represented in Figure 4.  

#### Table 3. Top 10 Significant Differentially Expressed Transcripts  
rank|target_id|gene name  
---|---|---  
1|ENSMUST00000075453.8|*Rpl21*  
2|ENSMUST00000002708.4|*Shh*  
3|ENSMUST00000031131.10|*Uchl1*  
4|ENSMUST00000031249.7|*Sparcl1*  
5|ENSMUST00000040576.9|*Parm1*  
6|ENSMUST00000056355.8|*Nat8l*  
7|ENSMUST00000058045.4|*Garem2*  
8|ENSMUST00000079324.13|*Ubl3*  
9|ENSMUST00000102553.10|*Hmgn2*  
10|ENSMUST00000112707.2|*Lrrc8b*  

#### Figure 4. Top 5 IGV Examples  
### Rp121  
![Rpl21](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Rpl21.PNG)  
### Shh  
![Shh](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Shh.PNG)  
### Uchl1  
![Uchl1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Uchl1.PNG)  
### Sparcl1  
![Sparcl1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Sparcl1.PNG)  
### Parm1  
![Parm1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Parm1.PNG)  

### Enhancer analysis  

First, we analyze the general signal of the histone modifications and the PhyloP track. H3K27ac peaks at the promoter of Rbm33, quickly falls off, and does not peak in Shh or other near gene regions. H3K4me1 peaks at the promoter of each gene in the region of interest, and gradually falls off, showing higher expression, particularly in intronic regions, compared to H3K27ac. The signal for these two histone modifications look largely similar between tissues, although we see some higher expression of MB-H3K4me1. This could be due to certain regulatory regions being specific to the mid-brain tissue, such that there is higher methylation to induce expression. Exonic regions seem to correlate directly with peaks in Phylop scores, which follows biological sense, as these regions are necessary to create protein product important for the organism. There are some high Phylop scores in intronic regions, and I can hypothesize that these are pivotal enhancer, activator, or otherwise critical regulatory regions without which the gene cannot create enough protein product to complete its function.  

Next, we look specifically at the ZRS region, which we know is implicated in regulation of *Shh* and potentially limb development. I visualized the results in Figure 5 and quantified them in Table 6. For FL and HL, the methylation levels were consistently high and the acetylation levels were consistently low, whereas for ML both methylation and acetylation were consistently low. The Phylop score is consistently high throughout this region, indicating that it is highly conserved.  

#### Figure 5. ZRS Region RNA-Seq, Histone Modifications, and PhyloP Scores  
![region of interest](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/regionOfInterest.PNG)  

#### Table 6. Quantification of Modification Scores  
Modification | Score Range | Scale  
---|---|---  
FL-H3K4me1|1.2 - 1.5 |0 - 3.5  
FL-H3K2ac|0.4 - 0.6|0 - 7.5  
HL-H3K4me1|1.3 - 1.6|0 - 4.4  
HL-H3K2ac|1.5 - 3.5|0 - 14.5  
ML-H3K4me1|0.3 - 0.5|0 - 3.5  
ML-H3K2ac|0.3 - 0.5|0 - 12.4  

Finally, we can analyze the MSA results from the outputted HTML file. I found two regions where there are missing sequences for snakes and conserved present sequences for all other organisms analyzed. These two regions are represented in Figure 7, where the last four sequences are the four snake organisms analyzed, and the other sequences are non-snakes.  

#### Figure 7. MSA Visualization of Absent Sequences in Snakes Exclusively  
![1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/zrs_1.PNG)  
![2](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/zrs_2.PNG)  

### Extra Credit  

1. Visualize expression of differentially expressed genes as a heatmap (e.g. with transcripts as rows and samples as columns). Cluster the rows by row and column. Do replicates cluster together? Are there clear clusters of up vs. down regulated genes in each tissue?  

To get all the TPM values, such that each sample is a column and each row is a transcript, we run a [custom script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/custom_scripts) (line 5). Then, we run a [python script](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/scripts/ec_1.py) to create a clustermap, which clusters the matrix by row and column. This cluster map is visualized in Figure 8 below. FL_Rep1 and FL_Rep2 tend to cluster, and MB_Rep1 and MB_Rep2 tend to cluster, whereas HL_Rep2 tends more towards clustering with MB compared to HL_Rep1. There are clear clusters of up vs. down regulated genes in each tissue, as in the clustermap below red represents up-regulation and blue represents down-regulation. It is evident that the red and the blue have separated into clear clusters.  

#### Figure 8. Clustermap of Transcripts and Samples  
![clustermap](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/tpm_clustermap.png)  

2. Perform gene ontology (GO) analysis or gene-set enrichment analysis (GSEA) on the top set of genes. There are many online tools for doing this (e.g. DAVID, Panther). What types of biological processes are enriched in differentially expressed genes? Inlude a table in your results section.  

Using [**DAVID Functional Annotation Tool (v6.8)**](https://david.ncifcrf.gov/summary.jsp)<sup>[**10**]</sup>, I performed gene-set enrichment analysis on the top 10 differentially-expressed genes and received the results in Table 9. The biological processes enriched in these genes involve structure development and cell proliferation, as well as surface proteins such as lipoproteins, glycoproteins, and membrane proteins.  

#### Table 9. Functional Annotation of Top 10 Differentially Expressed Genes  
![funcannotation](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Functional_Annotation_Chart.PNG)  

## Discussion  

In the ZRS region, H3K4Me1 levels are high for FL and HL, and from lecture slides as well as [studies on this modification](https://www.ncbi.nlm.nih.gov/pubmed/29255264)<sup>[**5**]</sup>, we know H3K4me1 plays an active role at enhancer regions by easing binding of chromatin regulators to control gene expression. Because we see high levels for only limb tissue and not brain tssue, we can hypothesize that the ZRS region is a putative enhancer region for the forelimb and hindlimb tissue, and therefore it is likely a critical region involved in limb development.  

There are other enhancer regions that seem to be limb-specific and are highly conserved. One example of this is chr5:29,212,459-29,212,837, where we see a high PhyloP score, high methylation for limb tissues, and low histone modifications for the midbrain tissue. To identify these regions computationally, we first set a range of the sequence to explore that qualifies as "nearby *Shh*". Then, we can computationally set thresholds for what constitutes a "high Phylop region", i.e. a region of the genome such that the PhyloP score is greater than 4.0 for a certain amount of bp (to indicate consistency). We can then subsect our region to just these highly conserved parts of the genome, and similarly computationally search for where in these regions we see a high (> some threshold) histone modification score for the limb tissues and a low (< some threshold) histone modification score for midbrain tissue.  

Alternatively, we can have a bar graph for these highly conserved regions of the genome, with the histone modification scores for midbrain tissues subtracted from the sum of the histone modification scores from the limb tissues as the y-axis, and then look at the peaks in the graph (> some threshold).  

Because the ZRS mutations we discovered are in the intronic region of a gene nearby *Shh*, we know that the mutations cannot necessarily mutate the amino acids of the protein product. Hence, what is more likely is that the mutations interfere with H3K4me1's ability to facilitate binding of chromatin regulators (such as the BAF complex) to the enhancer region. Maybe the mutation in snakes deletes a critical transcription factor binding site, or otherwise affects the ability for the enhancer to fully express *Shh*. As such, *Shh* cannot complete its role in early development of limb tissue, and snakes do not have limbs. In order to confirm these hypotheses, we would need to delete the sequences found to be absent in only snakes in the embryo of an animal with limbs, and monitor histone modifications of the ZRS, transcript and protein production levels of *Shh*, and phenotype of the animal as it develops.  

## Citations  

1. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
2. James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011)
3. Katoh, Misawa, Kuma, Miyata 2002 (Nucleic Acids Res. 30:3059-3066) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. (describes the FFT-NS-1, FFT-NS-2 and FFT-NS-i strategies)
4. Lettice, Laura A. et al. “The Conserved Sonic Hedgehog Limb Enhancer Consists of Discrete Functional Elements That Regulate Precise Spatial Expression.” Cell Reports 20.6 (2017): 1396–1408. PMC. Web. 2 May 2018.  
5. Local, A., Huang, H., Albuquerque, C.P., Singh, N., Lee, A.Y., Wang, W., Wang, C., Hsia, J.E., Shiau, A.K., Ge, K., Corbett, K.D., Wang, D., Zhou, H. & Ren, B. (2018) Identification of H3K4me1-associated proteins at mammalian enhancers. Nat. Genet. 50(1):73-82.
6. Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519
7. N P Brown, C Leroy, C Sander; MView: a web-compatible database search or multiple alignment viewer., Bioinformatics, Volume 14, Issue 4, 1 January 1998, Pages 380–381, https://doi.org/10.1093/bioinformatics/14.4.380
8. P.J. Kersey, J.E. Allen, A. Allot, M. Barba, S. Boddu, B.J. Bolt, D. Carvalho-Silva, M. Christensen, P. Davis, C. Grabmueller, N. Kumar, Z. Liu, T. Maurel, B. Moore, M. D. McDowall, U. Maheswari, G. Naamati, V. Newman, C.K. Ong, D.M. Bolser., N. De Silva, K.L. Howe, N. Langridge, G. Maslen, D.M. Staines, A. Yates. Ensembl Genomes 2018: an integrated omics infrastructure for non-vertebrate species Nucleic Acids Research 2018 46(D1) D802–D808, https://doi.org/10.1093/nar/gkx1011
9. UCSC Table Browser: Karolchik D, Hinrichs AS, Furey TS, Roskin KM, Sugnet CW, Haussler D, Kent WJ. The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 2004 Jan 1;32(Database issue):D493-6.
10. Glynn Dennis Jr., Brad T. Sherman, Douglas A. Hosack, Jun Yang, Michael W. Baseler, H. Clifford Lane, Richard A. Lempicki.  DAVID: Database for Annotation, Visualization, and Integrated Discovery. Genome Biology 2003 4(5): P3.  
11. Harold J. Pimentel, Nicolas Bray, Suzette Puente, Páll Melsted and Lior Pachter, Differential analysis of RNA-Seq incorporating quantification uncertainty, Nature Methods (2017), advanced access http://dx.doi.org/10.1038/nmeth.4324.
