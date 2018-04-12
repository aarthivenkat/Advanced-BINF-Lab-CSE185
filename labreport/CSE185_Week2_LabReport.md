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

#### Alignment of roommate data to reference FASTA  
We [download the reference FASTA file from NCBI](https://github.com/cse185-sp18/cse185-week2-aarthivenkat/blob/master/labreport/efetch_reference.txt).  

Now, using **[BWA-MEM](http://bio-bwa.sourceforge.net/) (v.0.7.12-r1039): Burrows-Wheeler Aligner**, we indexed the FASTA file with bwa index and ran the algorithm for alignment with bwa mem on the FASTA file and the two trimmed FASTQ files. This was outputted to a SAM file, so we can use **[SAMtools](http://samtools.sourceforge.net/) (v.1.5)** on the output SAM file to interpret percentage of reads aligned, create a BAM file. Using the samtools view -f 4 command, determine that 3430 reads are unmapped, so 283309 reads were mapped. We then make a pileup file using the FASTA file and BAM file in order to use VarScan to predict variants and effects.  

#### Common variants from VarScan  

We can now run **[VarScan](http://varscan.sourceforge.net) (v.2.3.9)** on the pileup file with a threshold of 95% and a --variants flag and --output-vcf flag to create a VCF file with the variant information. Finally, we would like to use a variant effect predictor to understand the potential implications of these variants, so we use an **[awk script](https://github.com/cse185-sp18/cse185-week1-aarthivenkat/tree/master/labreport/awk_script.txt)** to change the first column into the required format for the web-tool, and then use the **[Ensembl Variant Effect Predictor Tool](http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Tools/VEP)** for E.coli to determine whether the variant is in a gene (if so, which one), and if the mutation is missense, synonymous, or nonsense.  

This section should contain sufficient information so that other bioinformaticists could reproduce your results. You should briefly describe your raw data (what is it, what is the name of the reference) and describe what you did with it. You should write this in 2-3 paragraphs, not in a list. When you use a bioinformatics software program, do not write out the full command you typed, but do specify which program (ie ‘bwa-mem’ or ‘samtools tview’) you used and whether you used the default options. If you did not use the defaults, you should specify the exact settings you used. The first time you mention bioinformatics software or an online tool, you should cite it and specify which version of the tool you used. The correct citation for most software can be found by looking up its documentation on line (you don’t have to cite common tools like python or perl or the bash shell). If you write a custom script, (for example, our awk script from week 1), include that code in the labreport folder and reference it in your writeup.

## Results
This section should include the results of your data processing and data analysis, and may include tables with read lengths, pictures of quality distributions, or tables of gene names for examples. In the text, briefly restate how you got the results in full sentences, but in less detail than the methods, before you say what the results are (ie ‘reads were mapped to the reference and scanned to identify positions that likely contained mutations. We found….’). Refer to tables and figures by number, and include a brief descriptive title for each. Be sure to include any results specifically requested in the lab project tutorial. The results section should be as objective as possible, so please refrain from interpreting the meaning or significance here. It should be just the facts.

## Discussion
In 2-3 paragraphs, explain what you think the results mean, and why you are interpreting them this way. If you encountered any problems, or answered questions, discuss them and suggest ways to solve them with future experiments or analyses. Also include any information specifically requested in the tutorial.

## Citations
You can use any commonly used format you like, but be consistent. Lab reports will be submitted via turnitin to check for plagiarism, so be sure to cite other people’s ideas, and put everything in your own words (paraphrasing) if you aren’t using direct quotes.
