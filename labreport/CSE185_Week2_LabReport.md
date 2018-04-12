# Put title here
##### your name

## Abstract
In no more than 100 words, briefly summarize what was done in the lab this week, what the findings were, and why they were important.

## Introduction
In 2-3 paragraphs, provide enough background information to understand the biology behind the weeks project. Be sure to state what problem or question the week’s lab work addressed, and why it is important. You must cite at least one scientific journal article for this section (it can, but doesn’t have to be, the assigned reading). When you use outside resources, use in-text citations in the text attributing any ideas or information from materials outside of our course lecture or tutorial. In-text citations give the source for information right where it is written (1).
 
## Methods  

### Introducing Sequencing Data  

#### roommate.fastq  
I was given my roommate’s viral FASTQ file, which I was told to be the output of Illumina Single-end sequencing. There are **286739 reads** in this file, and each read was a different length, indicating that the data was pre-processed. This was confirmed with the medical school.  

To determine the number of cycles, we run an [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_roommate_cycles.txt) to output sorted unique read lengths (col2) and the count of each read length (col1).  

So, the maximum read length is 151 bp, so unless every read was trimmed, it is likely that the number of cycles carried out during the sequencing run is 151.  

#### reference fastq
To determine which of our "rare variants" are not solely sequencing errors, I use 3 control samples.  
SRR1705858.fastq 256586 reads  
SRR1705859.fastq 233327 reads  
SRR1705860.fastq 249964 reads  

The number of cycles during the sequencing run is 151, and the data was pre-processed because not every read has exactly 151 bp.   

### Alignment of roommate data to reference FASTA  
We [download the reference KF848938.1 file from NCBI in the format FASTA](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/efetch_reference.txt).  

Now, using **[BWA-MEM](http://bio-bwa.sourceforge.net/) (v.0.7.12-r1039): Burrows-Wheeler Aligner**, we indexed the FASTA file with bwa index and ran the algorithm for alignment with bwa mem on the FASTA file and the two trimmed FASTQ files. This was outputted to a SAM file, so we can use **[SAMtools](http://samtools.sourceforge.net/) (v.1.5)** to create a BAM file. Using the samtools view -f 4 command, determine that 3430 reads are unmapped, so 283309 reads were mapped. We then make a pileup file using the FASTA file and BAM file in order to use VarScan to predict variants and effects.  

### Common variants from VarScan  

We can now run **[VarScan](http://varscan.sourceforge.net) (v.2.3.9)** on the pileup file with a threshold of 95% and a --variants flag and --output-vcf flag to create a VCF file with the variant information. To format the VCF file for readability, we run an [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_vcf_interpret.txt) to output common variants. 

Finally, we would like to manually check to see how the amino acids change to determine potential implications of these variants. We copy the reference FASTA sequence into online sequence editor [WebDSV](http://www.molbiotools.com/WebDSV/) and note which variants are missense or synonymous, as well as the amino acid change information.  

### Rare variants from VarScan

Similarly, we run VarScan on the pileup with a threshold of 0.1% and a --variatns and --output-vcf flag to create a VCF file with the rare variant information. We want to format the VCF with frequency this time, so run a different [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_vcf_interpret_freq.txt) to output the initial list of "variants" with frequency above 0.1%.  

### Sequencing Errors in Control Files  

To determine the sequencing errors in the control files, we use the exact same process as above to find the rare variants. We align the control fastq files to the reference and create the mpileup file for VarScan using BWA MEM and SAMTools. Then we run VarScan with minimum variant frequency 0.1%, only outputting variants and formatting in VCF. Then, I parse the VCF using the [awk script](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/awk_vcf_interpret_freq.txt) from the rare variants section above.  

### Excel to Compare Roommate to Control  

We create an [Excel spreadsheet](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/control_sample_variants.xlsx) with all the errors of the control files as well as the errors and variants of the roommate file. Using Excel commands `AVERAGE` and `STDEV`, we get the average and stdev for each of the control files, and then use a [nested if statement](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/excel_greaterthan3std.txt) for each of the roommate variants to determine if the frequency is greater than 3 standard deviations from the average of any one of the control files.  

We consider this list of roommate variants to be "real variants", as the frequency is significantly greater than the frequency of sequencing errors. We plug these variants again into WebDSV and note which variants are missense or synonymous, as well as the amino acid change information, which is critical for epitope analysis.  

### Epitope analysis  

Simply reading the [epitope paper given to us](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4482133/#R19), I manually scanned the results section to see which residues in the epitopes were implicated with antigenic variation, and whether any of these residues coincided with residues from the roommate's rare variants.  

## Results
This section should include the results of your data processing and data analysis, and may include tables with read lengths, pictures of quality distributions, or tables of gene names for examples. In the text, briefly restate how you got the results in full sentences, but in less detail than the methods, before you say what the results are (ie ‘reads were mapped to the reference and scanned to identify positions that likely contained mutations. We found….’). Refer to tables and figures by number, and include a brief descriptive title for each. Be sure to include any results specifically requested in the lab project tutorial. The results section should be as objective as possible, so please refrain from interpreting the meaning or significance here. It should be just the facts.

## Discussion
In 2-3 paragraphs, explain what you think the results mean, and why you are interpreting them this way. If you encountered any problems, or answered questions, discuss them and suggest ways to solve them with future experiments or analyses. Also include any information specifically requested in the tutorial.

## Citations
You can use any commonly used format you like, but be consistent. Lab reports will be submitted via turnitin to check for plagiarism, so be sure to cite other people’s ideas, and put everything in your own words (paraphrasing) if you aren’t using direct quotes.
