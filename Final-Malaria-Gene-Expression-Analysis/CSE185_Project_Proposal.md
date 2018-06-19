# Transcriptome Expression Analysis of Chloroquine-Resistant Plasmodium Falciparum
#### Aarthi Venkat
#### 24 May 2018

## Biological Question

The emergence of drug resistance in malaria is a critical public health concern, and, for decades, researchers have employed variousc approaches to studying the *Plasmodium* parasitic species. I am interested in adding to the literature of the molecular mechanism of malaria drug resistance through a whole transcriptome analysis of *Plasmodium falciparum*. I will first replicate the processing pipeline from the raw sequencing reads to readable, differential expression output. Then, I plan to extend the statistical analysis of the researchers by utilizing *cummeRbund*, a tool for visualization and statistics of Cuffdiff output.  

## Dataset

The data is paired-end RNA Seq data from Illumina HiSeq 2500 sequencing. It is annotated with TopHat from the Ensemble Genome Database (Release 26). Reads were collected from three strains of *Plasmodium falciparum* - one non-resistant strain (Pfal3D7;SRX1546678) and two resistant, one cultured in the presence of chloroquine (PfalDd2_CQ;SRX1535489) and one in the absence (PfalDd2;SRX1546677). The raw data can be found in NCBI BioProject under the accession number **PRJNA308455**, and the processed data (Cuffdiff output) can be found in the GEO under the accession number **GSE77499**. Ideally, my pipeline will output expression data similar to the published processed output, such that downstream analysis renders reliable results. Each of the three raw files are in FASTQ format, and are each 1.8 GB, so I will perform differential analysis on some subset of the files. The processed file is in EXP.DIFF format, which is the standard output of Cuffdiff, expression estimation and differential expression tool. It is 650 KB in size.  

The publication linked to the data is ["Whole transcriptome expression analysis and comparison of two different strains of Plasmodium falciparum using RNA-Seq."](https://www.ncbi.nlm.nih.gov/pubmed/27222812)  

## Bioinformatics pipeline  

The paper's method of processing the RNA-Seq data was as follows:  

1. QC check with FastQC  
2. Base trimming through in-house Perl script  
3. Adapter trimming with Cutadapt  
4. Contamination removal with bowtie2, various in-house scripts, and picardtools  
5. Read alignment with Tophat and default parameters  
6. Expression estimation with Cufflinks  
7. Differential expression analysis with Cuffdiff  

I will replicate this pipeline in my own analysis, yet I will use Sickle to trim the bases instead of an in-house script, as I am comfortable with the tool and understand its usage. Furthermore, I will likely not include contamination removal, as the paper did not go into detail about which subtools of bowtie2 and picardtools the researchers used to remove contamination.  

The second step of analysis will involve expanding the analysis from the basic volcano plots the researchers created. They only briefly attempted to characterize differential expression, so I hope to incorporate cummeRbund, which I have already installed locally, to work with my cuffdiff output and perform comparative analyses such as PCA, gene ontology analysis, and identification of significant gene sets.  
