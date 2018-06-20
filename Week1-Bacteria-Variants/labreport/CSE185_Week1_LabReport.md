# What Causes Antibiotic Resistance? 
## Using Bioinformatics Sequencing Tools to Identify Key Variants in E.coli Ampicillin-resistant Genome
##### Aarthi Venkat

## Abstract  
This week, we leveraged research and industry-standard bioinformatics tools, including FastQC, Sickle, BWA-MEM, Samtools, and VarScan, to identify variants from raw data in an ampicillin-resistant strain of Escherichia coli. Through our analysis, we were able to implicate six genes (b0084, b0771, b1821, b3404, b4161), the first five corresponding to missense mutations and the last a synonymous mutation. We look to mechanisms of antibiotic resistance particular to these genes, and posit potential treatment recommendations. This paper hopes to stress the necessity of microbial sequencing as a diagnostic standard to healthcare.    

## Introduction
Antibiotics have known cellular targets and are designed to kill or impair bacteria-specific features, such as the bacterial cell wall,  cell membrane, essential enzymes, and protein synthesis pathways. As such, bacteria have developed antibiotic resistance to avoid antibiotics — developments as a result from mutations selectively advantaging certain bacteria for survival. Therefore, it is of interest to determine which genes and antibiotic resistance mechanisms may be present in a patient to ensure the highest possible accuracy for the lowest time and monetary cost.  

According to [Peacock](https://www.nature.com/news/health-care-bring-microbial-sequencing-to-hospitals-1.15282), with traditional methods of determining genes implicated in resistance, such as PCR, researchers can detect only a few resistance markers in a sample, and getting results can take weeks while the patient continues to take ineffective treatment. Microbial sequencing and variant analysis can not only be inordinately faster and contain a broader reach of results and statistics, but also be adopted into clinical diagnostics to determine potential treatments for the patient instantaneously.  

Thus, this paper looks to leverage sequencing software and variant analysis tools to determine potential genes implicated in ampicillin resistance, so as to understand treatment targets for this strain through bioinformatics approaches exclusively. Our reference genome will be Escherichia coli str. K-12 substr. MG1655, complete genome.  

_________________________________________________________________________________________
## Methods
After studying the given FASTQ files (amp_res_1.fastq and amp_res_2.fastq) manually and acquainting myself with the raw data, I ran **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v.0.11.7): High-Throughput Sequence QC Analysis Tool** on the two FASTQ files (from paired-end sequencing). This produced HTML output, which allowed me to see the failures of the raw data and opt for trimming of the ends to improve per base sequence quality. So, I turned to **[Sickle](https://github.com/najoshi/sickle) (v.1.33): Windowed Adaptive Trimming for FASTQ files** and trimmed the reads, cutting off where sequence quality was below the set threshold. Using sickle pe (for paired end) and sanger as the quality type, we tried two thresholds — the default, 20, and 30 — and used the default threshold in subsequent analyses.  

Now, using **[BWA-MEM](http://bio-bwa.sourceforge.net/) (v.0.7.12-r1039): Burrows-Wheeler Aligner**, we indexed the FASTA file with bwa index and ran the algorithm for alignment with bwa mem on the FASTA file and the two trimmed FASTQ files. This was outputted to a SAM file, so we can use **[SAMtools](http://samtools.sourceforge.net/) (v.1.5)** on the output SAM file to interpret percentage of reads aligned, create a BAM file, and tview the BAM file and FASTA file to visualize the alignments interactively. We then make a pileup file using the FASTA file and BAM file in order to use VarScan to predict variants and effects.  

Downloading **[VarScan](http://varscan.sourceforge.net) (v.2.3.9)** from sourceforge, we run the program on the pileup file with a threshold of 70% to create a VCF file with the variant information. Finally, we would like to use a variant effect predictor to understand the potential implications of these variants, so we use an **[awk script](https://github.com/cse185-sp18/cse185-week1-aarthivenkat/tree/master/labreport/awk_script.txt)** to change the first column into the required format for the web-tool, and then use the **[Ensembl Variant Effect Predictor Tool](http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Tools/VEP)** for E.coli to determine whether the variant is in a gene (if so, which one), and if the mutation is missense, synonymous, or nonsense.  

_________________________________________________________________________________________  
## Results

In order to understand how our data was manipulated through quality control and alignment, we can track the number of reads in Table 1.  

### Table 1. Number of Reads  

Initial number of reads | Number of reads after trimming | Number aligned
---|---|---
7107040|6904494|6899109 (primary and alternative mappings included)  

amp_res_1.fastq and amp_res_2.fastq have 3553520 reads in each file, meaning there were initially 7107040 reads. After quality trimming using Sickle, there were 6904494 reads left. Using samtools flagstat, we know that 6899109 of those reads were aligned (primary and alternative mappings included).  

Using FastQC, we were able to visualize the per base sequence quality in HTML. Over the run, the quality degrades, so we run Sickle to trim the poorly-sequenced ends of the sequence and rerun FastQC. In Figure 2, we include the visualizations of the sequence before and after trimming.  

### Figure 2. Per Bases Sequence Quality Before and After Trimming  

#### amp_res_1_before_trim  
![](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/tree/master/Week1-Bacteria-Variants/labreport/res1_beforeTrim.PNG)
#### amp_res_2_before_trim  
![](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/tree/master/Week1-Bacteria-Variants/labreport/res2_beforeTrim.PNG)
#### amp_res_1_after_trim  
![](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/tree/master/Week1-Bacteria-Variants/labreport/res1_afterTrim.PNG)
#### amp_res_2_after_trim  
![](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/tree/master/Week1-Bacteria-Variants/labreport/res2_afterTrim.PNG)

Evidently, trimming the ends of the sequences was able to improve per base sequence quality, allowing us to move forward with subsequent analyses.  

After trimmming, we ran BWA-MEM to align the reads to the reference genome, and then SAMtools and create a pileup file to interpret the alignment. Then, we can use VarScan on the pileup file to extract interesting information. We find that six positions are above the minimum threshold of variant allele frequency set (70%).  

VarScan.vcf has the following information indicating mutations in the E.coli strain. I have recorded the position, reference, alternative base, and variant allele frequency in each column in Table 3 below.  

### Table 3. Variant Base and Frequency Information  

POS | REF BASE | ALT BASE | FREQ
---|---|---|---
92439|G|A|99.32%
803662|C|A|100%
852762|A|G|100%
1905761|A|G|100%
3535147|G|A|100%
4390754|G|T|100%

Using the Variant Effect Predictor tool for E.coli, and returning results for variants in coding regions only, we begin to implicate five genes, listed in Table 4.  

### Table 4. Variant Effect Predictions and Genes Implicated  

POS | OCCURS IN GENE? | TYPE | GENE NAME | AA SUBSTITUTION
---|---|---|---|---
92439|yes|missense|b0084|D to N
803662|yes|missense|b0771|L to I
852762|no|N/A|N/A|N/A
1905761|yes|missense|b1821|G to D
3535147|yes|missense|b3404|V to G
4390754|yes|synonymous|b4161|A to A  

Four genes are involved in missense mutations. Each gene has a particular metabolic function that could be implicated in antibiotic resistance. Table 5 is the information on each gene collated from [EcoCyc](https://ecocyc.org/).  

### Table 5. Gene and Protein Function  

Gene | Protein | Protein function | Potential antibiotic mechanism of mutation
---|---|---|---
b0084|Ftsl|essential protein for bacterial cell-division and a transpeptidase|resists interference of antibiotic with cell wall creation
b0771|YbhJ|oxidative stress response|mutate to avoid antibiotic-induced oxidative stress
b1821|mntP protein|manganese export protein|mutate to fortify efflux of manganese component of an antibiotic
b3404|EnvZ|kinase regulating expression of porin genes| mutation may decrease permeability of the porin genes to block drug from entering the cell  

The potential mechanisms of resistance are further involved in the discussion.  

_________________________________________________________________________________________
## Discussion  

### Four Genes of Interest  

After interpreting the data from FASTQ to VCF file, we have a set of genes of which we could biologically study the function and mutation. Before this, it is necessary to understand the mechanisms of antibiotic resistance common to bacterial pathogens, and how these genes could potentially be involved in resistance. According to lecture 2, the four mechanisms of antibiotic resistance are as follows:

1. Target site alteration, prohibiting target binding  
2. Inactivation or modification of the antibiotic to render it ineffective  
3. Alteration of the pathway to detour the site of antibiotic attack  
4. Reduction of the amount of the drug in the cell through either decreased permeability of the membrane or efflux pumps (effectively kicking the drug out)  

As such, we can make predictions as to what the mutations in three of the genes may be doing.  

[**b0084**](https://biocyc.org/gene?orgid=ECOLI&id=EG10341) codes for FtsI, an essential protein for bacterial cell-division and a transpeptidase, according to BioCyc. Thus, as discussed in lecture, its function is to convert peptidoglycan into cross-linked peptidoglycan to create the bacterial cell wall. Ampicillin, the antibiotic of interest in this lab, impairs bacteria by binding to Ftsl, interfering with the cell wall creation pathway and rendering the cell wall broken. Thus, a mutation in the gene b0084 may be leveraging the first mechanism of antibiotic resistance — target site alteration — to create a working protein Ftsl with a target site alteration, so that the antibiotic cannot bind to Ftsl and Ftsl can complete its function.  

[**b1821**](https://biocyc.org/gene?orgid=ECOLI&id=G6999) codes for a manganese export protein, as per BioCyc. According to [Fisher et al](http://jb.asm.org/content/198/20/2810.full?cited-by=yes&legid=jb;198/20/2810#cited-by), manganese plays a critical role in cellular physiology and metabolism of bacteria, but the amount of manganese in the cell must be highly regulated, as the accumulation of toxic levels of manganese can disrupt cellular function. If the antibiotic of interest has manganese as a major component, such as the manganese tricarbonyl complex discussed by [Betts et al](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0186359#), to the point where there is a major influx of manganese in the cell, a mutation in the manganese export protein which fortifies its function could potentially efflux the drug out of the cell, leveraging the fourth mechanism in the list above.  

[**b3404**](https://biocyc.org/gene?orgid=ECOLI&id=EG10269) codes for a kinase with the function of regulating expression of the outer membrane porin genes (ompF and ompC) in response to changes in osmotic pressure, according to BioCyc. As discussed in lecture, a mutation in the kinase gene could allow it to increase regulation of expression of the porin genes, to the point of decreasing permeability of the membrane and disallowing the drug to enter the cell.  

### Treatment Recommendation  

As discussed in lecture, if we implicate b0084, which uses the mechanism of target site alteration, then it would be in the patient’s best interest to try a secondary treatment with a new target. If b1821 is playing a significant role in antibiotic resistance, then we must target the efflux pumps by prescribing antibiotics with manganese efflux pump inhibitors. Finally, if we need to target b3404, we must prescribe antibiotics with permeabilizers to ensure the drug makes it into the cell.  

### Problems and Questions Encountered  

In the first round of quality control with FastQC, we discovered errors for Per Base Sequence Quality, Per Tile Sequence Quality, and Per Base Sequence Content. While we reconciled the error of per base sequence quality by trimming, it would be worth investigating why the other two errors occurred, and whether these errors would be detrimental to subsequent analyses.  

Further, there were several places where the user was required to place somewhat arbitrary thresholds; for example, the quality threshold for sickle was set to 20, and the variant allele frequency threshold set to 70%. It would be wise to perform the same analysis in the future with multiple thresholds and look for significant differences.  

### What do you think you learned?  

Through whole genome sequencing analysis, I was able to — very quickly — determine potential genes of interest, their functions, and potential implications of how they developed antibiotic resistance. This would not have been possible using traditional non-computational methods, which can detect only a handful of resistance markers in a sample. With the introduction of sequencing into local, clinical diagnostics, as well as the creation of automated tools to interpret these genomes, we can be able to control infections — as long as we push to leverage these computational advances.  

## Citations  
1. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc  
2. Betts J, Nagel C, Schatzschneider U, Poole R, La Ragione RM (2017) Antimicrobial activity of carbon monoxide-releasing molecule [Mn(CO)3(tpa-κ3N)]Br versus multidrug-resistant isolates of Avian Pathogenic Escherichia coli and its synergy with colistin. PLoS ONE 12(10): e0186359. https://doi.org/10.1371/journal.pone.0186359  
3. Fisher CR, Wyckoff EE, Peng ED, Payne SM. 2016. Identification and characterization of a putative manganese export protein in Vibrio cholerae. J Bacteriol 198:2810–2817. doi:10.1128/JB.00215-16.
4. Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files 
(Version 1.33) [Software].  Available at https://github.com/najoshi/sickle.  
5. Keseler, I.M., Mackie, A., Santos-Zavaleta, A., Billington, R., Bonavides-Martinez, C., Caspi, R., Fulcher, C., Gama-Castro, S., Kothari, A., Krummenacker, M., Latendresse, M., Muñiz-Rascado, L., Ong, Q., Paley, S., Peralta-Gil, M., Subhraveti, P., Velazquez-Ramirez, D.A., Weaver, D., Collado-Vides, J., Paulsen, I., and Karp, P.D., The EcoCyc database: reflecting new knowledge about Escherichia coli K-12, Nucleic Acids Research 45:D543-550 2017.  
6. Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111  
7. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]  
8. Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]  
9. Peacock, Sharon. “Health Care: Bring Microbial Sequencing to Hospitals.” Nature News, Nature Publishing Group, 22 May 2014, www.nature.com/news/health-care-bring-microbial-sequencing-to-hospitals-1.15282.  
10. P.J. Kersey, J.E. Allen, A. Allot, M. Barba, S. Boddu, B.J. Bolt, D. Carvalho-Silva, M. Christensen, P. Davis, C. Grabmueller, N. Kumar, Z. Liu, T. Maurel, B. Moore, M. D. McDowall, U. Maheswari, G. Naamati, V. Newman, C.K. Ong, D.M. Bolser., N. De Silva, K.L. Howe, N. Langridge, G. Maslen, D.M. Staines, A. Yates. Ensembl Genomes 2018: an integrated omics infrastructure for non-vertebrate species
Nucleic Acids Research 2018 46(D1) D802–D808, https://doi.org/10.1093/nar/gkx1011  
