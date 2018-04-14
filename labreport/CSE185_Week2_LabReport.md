# Put title here
##### your name

## Abstract
In no more than 100 words, briefly summarize what was done in the lab this week, what the findings were, and why they were important.

## Introduction
In 2-3 paragraphs, provide enough background information to understand the biology behind the weeks project. Be sure to state what problem or question the week’s lab work addressed, and why it is important. You must cite at least one scientific journal article for this section (it can, but doesn’t have to be, the assigned reading). When you use outside resources, use in-text citations in the text attributing any ideas or information from materials outside of our course lecture or tutorial. In-text citations give the source for information right where it is written (1).
 
## Methods  

### Introducing Sequencing Data  

#### roommate FASTQ file  
I was given my roommate’s viral FASTQ file, which I was told to be the output of Illumina Single-end sequencing. Each read in this file is a different length, indicating that the data was pre-processed. This was confirmed with the medical school.  

To determine the number of cycles, we ran an [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_roommate_cycles.txt) to output sorted unique read lengths (col2) and the count of each read length (col1).  

The maximum read length is 151 bp, so unless every read was trimmed, it is likely that the number of cycles carried out during the sequencing run is 151.  

#### reference FASTQ files  
To determine which of our "rare variants" are not solely sequencing errors, I use 3 control samples (SRR1705858, SRR1705859, SRR1705860).  

The number of cycles during the sequencing run is 151, and the data was pre-processed because not every read has exactly 151 bp.   

### Alignment of roommate data to reference FASTA  
We first [download the reference KF848938.1 file from NCBI in the format FASTA](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/efetch_reference.txt).  

Then, using **[BWA-MEM](http://bio-bwa.sourceforge.net/) (v.0.7.12-r1039): Burrows-Wheeler Aligner**, we index the FASTA file with bwa index and run the alignment with bwa mem on the FASTA file and the roommate FASTQ file. This was outputted to a SAM file, so we can use **[SAMtools](http://samtools.sourceforge.net/) (v.1.5)** to create a BAM file. We then make a mpileup file using the FASTA file and BAM file in order to use VarScan to predict variants and effects.  

### Common variants from VarScan  
We can now run **[VarScan](http://varscan.sourceforge.net) (v.2.3.9)** on the mpileup file with a threshold of 95% and a --variants flag and --output-vcf flag to create a VCF file of common variants. To format the VCF file for readability, we run an [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_vcf_interpret.txt) to output common variants. 

Finally, we would like to manually check the file to determine potential implications of these variants. We copy the reference FASTA sequence into online sequence editor [WebDSV](http://www.molbiotools.com/WebDSV/) and note which variants are missense or synonymous, as well as the amino acid change information.  

### Rare variants from VarScan  
Similarly, we run VarScan on the mpileup with a threshold of 0.1% and a --variants and --output-vcf flag to create a VCF file with the rare variant information. We want to format the VCF with frequency this time, so run a different [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_vcf_interpret_freq.txt) to output the initial list of "variants" with frequency above 0.1%.  

### Sequencing Errors in Control Files  
To determine the sequencing errors in the control files, we use the exact same process as above to find the rare variants. We align the control fastq files to the reference and create the mpileup file for VarScan using BWA MEM and SAMTools. Then we run VarScan with minimum variant frequency 0.1%, only outputting variants and formatting in VCF. Then, I parse the VCF using the [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_vcf_interpret_freq.txt) from the rare variants section above. Because there will be no variants in the control files, any "variants" we do find will be sequencing errors with frequency >0.1%.  

### Excel to Compare Roommate to Control  
We create an [Excel spreadsheet](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/all_sample_variants.xlsx) with all the errors of the control files as well as the errors and variants of the roommate file. Using Excel commands `AVERAGE` and `STDEV`, we get the average and stdev for each of the control files, and then use a simple [if statement](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/excel_greaterthan3std.txt) for each of the roommate variants to determine if the frequency is greater than 3 standard deviations from the average of any one of the control files.  

We consider this list of roommate variants to be "real variants", as the frequency is significantly greater than the frequency of sequencing errors. We plug these variants again into WebDSV and note which variants are missense or synonymous, as well as the amino acid change information, which is critical for epitope analysis.  

### Epitope Analysis  
Simply reading the [epitope paper given to us](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4482133/), I manually scanned the results section to see which residues in the epitopes were implicated with antigenic variation, and whether any of these residues coincided with residues from the roommate's rare variants.  

## Results  
This section should include the results of your data processing and data analysis, and may include tables with read lengths, pictures of quality distributions, or tables of gene names for examples. In the text, briefly restate how you got the results in full sentences, but in less detail than the methods, before you say what the results are (ie ‘reads were mapped to the reference and scanned to identify positions that likely contained mutations. We found….’). Refer to tables and figures by number, and include a brief descriptive title for each. Be sure to include any results specifically requested in the lab project tutorial. The results section should be as objective as possible, so please refrain from interpreting the meaning or significance here. It should be just the facts.  
First, we will look at the mappability for each of the four samples (Table 1).  

### Table 1. Reads Mapped
Sample|Initial Number of Reads| Number of Reads Mapped  
---|---|---  
roommate|286739|283309  
SRR1705858|256586|256500  
SRR1705859|233327|233251  
SRR1705860|249964|249888  

For each control file, we can determine the error rate by calculating the average/stdev of all the variant allele frequencies (Table 2).  

### Table 2. Control Sample Stats  
Control|AvgFreq|StdevFreq  
---|---|---  
SRR1705858|0.256%|0.0717%  
SRR1705859|0.237%|0.0524%  
SRR1705860|0.250%|0.0780%  

From VarScan and WebDSV, we were able to determine the implications of common variants for the roommate file. Similarly, after choosing the subset of roommate rare variants that were greater than 3 standard deviations from any one of the averages of the control files, we used WebDSV on that subset to add to the following table (Table 3).   

### Table 3. Roommate Variants  
Variant | Original Codon > Mutated Codon | AA Change | Type of Mutation | Type of Variant  
---|---|---|---|---  
A72G | ACA > ACG| Thr24Thr | synonymous | common  
C117T | GCC > GCT | Ala39Ala | synonymous | common  
G595T | GCA > TCA | Ala199Ser | missense | common  
T774C | TTT > TTC | Phe258Phe | synonymous | common  
T1008G | GCT > GCG | Ala336Ala | synonymous | common  
A1260C | CTA > CTC | Leu420Leu | synonymous | commmon  
T1339C | TTG > CTG | Leu447Leu | synonymous | common   
T38C | CTG>CCG | Leu13Pro | missense | rare  
C495T | AAC>AAT | Asn165Asn | synonymous | rare  
G910A | GCC>ACC| Ala304Thr | missense | rare  
G1293A |CTG>CTA | Leu431Leu | synonymous | rare  
G1521A |CTG>CTA | Leu507Leu | synonymous | rare  

### Epitope Analysis  

Finally, we turn to Munoz and Deem to determine whether these residues are located in epitopes of hemagglutinin. Of the residues in Table 3 (column 3), only **residues 165 and 304** were in epitopes B and C (respectively). Because reside 165 is synonymous and 304 is missense, it is likely that only 304 is implicated in this case, and **a mutation in epitope C is a factor in my getting the flu**.  

## Discussion
In 2-3 paragraphs, explain what you think the results mean, and why you are interpreting them this way. If you encountered any problems, or answered questions, discuss them and suggest ways to solve them with future experiments or analyses. Also include any information specifically requested in the tutorial.

## Citations
You can use any commonly used format you like, but be consistent. Lab reports will be submitted via turnitin to check for plagiarism, so be sure to cite other people’s ideas, and put everything in your own words (paraphrasing) if you aren’t using direct quotes.
