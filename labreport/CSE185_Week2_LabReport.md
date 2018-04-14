# Why did I get the Flu?  
### Rare Missense Mutation in Epitope C implicated in Roommate Strain of Influenza  
##### Aarthi Venkat  

## Abstract  

We determine the cause of my flu knowing that my roommate's virus matches the HI profile for 2017/2018 strain H3N2. Suspecting that a portion of the population could have mutated such that the antigen could not be recognized by my immune system, we sequence my roommate's viral sample. Ultimately, we find seven common (residues 24, 39, 199, 258, 336, 420 447) and five rare (residues 13, 165, 304, 431, 507) variants. Residue 304 is both missense and located in an epitope (Epitope C), indicating that it may be implicated in antigenic variation and the cause of my flu.  

## Introduction  

According to the [CDC](www.cdc.gov/flu/), flu vaccines develop protection against viral infection. The vaccine allows the immune system to produce antibodies to recognize surface of hemagglutinin, the protein facilitating viral binding to respiratory or red blood cells. The vaccine is designed based on research indicating the most commmon viruses in the upcoming season. However, mutations in epitope regions on the surface of HA mean that the immune system can no longer recognize the viral species. Such an occurence is called antigenic variation, and is a common viral mechanism to evade immune response. As such, influenze mutates at a rate of one mutation per genome per replication, and multiple quasispecies can exist in a single host. Therefore, it is not uncommon that a person who got the flu vaccine could also get the flu, as mutations allowing the virus to escape the vaccine are a regular aspect of influenza. 

As such, we are interested in studying my roommate's viral DNA to determine which mutations in the viral DNA may have allowed the  virus to escape the flu vaccine. However, because quasispecies tend to exist in a single host, we must have extremly high coverage, known as deep sequencing, to be able to predict rare variants that the HI assay may miss.

While increasing coverage is a smart way to unveil rare variants, it brings with it another issue. Next generation sequencing is not a perfect method, and errors may be introduced in every part of the pathway, from sampling of the DNA to bioinformatics analyses. Acute efforts must be made to mitigate these problems, and in this lab we error control by using replicates as our controls, as encouraged by [Robasky et al](https://www.nature.com/articles/nrg3655), by increasing coverage for rare variant analysis, and by trimming the reads based on per base quality. We also determine the variant frequency in control samples, indicating a base error rate for sequencing, and then subset our rare variants as those with frequency significantly higher than the error rate.
 
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
T38C | CTG > CCG | Leu13Pro | missense | rare  
C495T | AAC > AAT | Asn165Asn | synonymous | rare  
G910A | GCC > ACC| Ala304Thr | missense | rare  
G1293A |CTG > CTA | Leu431Leu | synonymous | rare  
G1521A |CTG > CTA | Leu507Leu | synonymous | rare  

### Epitope Analysis  

Finally, we turn to Munoz and Deem to determine whether these residues are located in epitopes of hemagglutinin. Of the residues in Table 3 (column 3), only **residues 165 and 304** were in epitopes (B and C, respectively). Because reside 165 is synonymous and 304 is missense, it is likely that only 304 is implicated in this case, and **a mutation in epitope C is a factor in my getting the flu**.  

## Discussion

### Interpreting Results  

In Table 1, we see that all our sampes have largely the same number of reads, as well as high mappability. Both of these factors will help mitigate errors in terms of bias. Table 2 indicates that all our controls have an error rate of roughly 0.25%, with standard deviation rate of roughly 0.07%. Because all these error rates are nearly the same, we interpret the roommate data by subsetting only variants with frequency significantly higher than one of these error rates. Table 3 shows our final list of all the roommate variants we believe are significant, but because the only regions of the virus that are specific to immune system response are the epitope regions, we would like to further subset this list to those mutations are involved in antigenic variation, which we do in the epitope anlaysis component of our results.  

### How did you decide which mutations were most likely to be real?  

For common variants, we rely on the threshold of 95% to argue that, at this threshold, we see all the roommate sequences with this variation, and thus they are likely real because the threshold is significantly high.  

For rare variants, we could not simply rely on the low threshold of 0.1% because this includes variants introduced through faults in deep sequencing (in other words, errors). Thus, we had to determine the error rate of the sequencing experiments by finding the variants in the control files, and then only choose roommate variants with frequency significantly higher than this error rate. We decided that mutations satisfying this condition were also most likely real.  

### Explain how you think you were able to get the flu from your roommate, even though you had received the flu vaccine.  

According to Munoz and Deem, epitopes are components of the antigen recognized by the immune system, yet antigenic variation occurring in these regions remains an effective mechanism for viruses to evade adaptive response of the host immune system. As such, it is resourceful to study mutations in epitope regions to determine how the viral DNA has changed to escape immune response. Of the residues in Table 3 column 3, which shows significant variants in the roommate DNA, only two residues (165 and 304) were in epitopes, (B and C, respectively). Because reside 165 is synonymous and 304 is missense, it is likely that only 304 is implicated in this case, and a mutation in epitope C is changing the antigen conformation such that the body cannot recognize this new antigen, and it cannot respond to it - hence I get the flu. It would interesting to see if the residue 304 mutation affects immune response to epitope B, as B is a dominant epitope.  

### Error Control for Deep Sequencing 

As discussed in lecture, sequencing errors can occur in every part of the process - sample collection, library preparation, sequencing and analysis. The following laboratory steps can be taken to minimize errors.  

1. Sample collection: contamination can introduce DNA into the sample that is different than the source of analysis; this can result in gross misanalysis of the data, particularly because it occurs upstream of the sequencing pathway. One way to control for this would be to ensure all people involved in the sampling are acutely aware of lab protocol, such as wearing the proper equipment and following sterilization and preparation procedures.  

2. Sample collection: amplification mutations can change how the cluster is read by the sequencing machine, especially if they are introduced in an early amplification cycle such that the DNA will be amplified with the mutation. A way to control this process is to minimize amplification using PCR-free library prep, which necessarily requires more input DNA.  

3. Informatics: Using quality scores from the FASTQ files, we are able to understand how the quality lowers as the read length increases. This is largely because errors that occur during sequencing tend to accumulate, resulting in less certainty as to the base represented by a cluster. To ameliorate this problem, we can trim the reads to ensure the per base quality scores satisfies a threshold for every read. This was done for us in this lab. We can also use the per tile quality scores to determine if there were significant issues with the flow cell, such as a bubble or a piece of dirt.  

Error control is important for accurately identifying and analyzing rare variants because rare variants will have extremely low frequency - that close to the sequencing error rate. If we have too many errors in our sequencing, then the error rate will be inordinately high, and we will not be able to determine whether a variant we have identified is a true variant or an error. However, by applying principles of error control we learned in lecture and from [Robasky et al](https://www.nature.com/articles/nrg3655), we can lower error rate and increase the significance of our findings.  

### Problems and Questions  

It seems that our control files have a relatively high sequencing error rate (0.25%), so I would be interested in looking into why these controls consistently differ from the reference sequence at certain positions. If I knew the cause, I could further investigate how to mitigate these effects and analyze rare variants with stronger confidence. Further, the epitope analysis paper by Munoz and Deem seems to admit there are lacuns in the quanitification of epitope dominance, particularly that a "precise determination" of dominance and an experimental measure of p_epitope and cross activity in proposed vaccine strains would be productive. This paper seems to be an initial analysis of epitopes, and with more quantification and measures of epitope variation for future experiments, we could make better predictions with our data.  

## Extra Credit  

For dataset roommate, I would calculate actual average coverage by first putting all the mapped reads into a text file, and then determining the number of bp in that file, and finally dividing that number by 1665, the number of bp in the reference file.  

To put all the mapped reads into a text file, I would run `samtools view -F 4 roommate.bam | cut -f 10 > roommateMappedReads.txt`, which filters the reads in the roommate.bam file to those that are mapped, then grabs the 10th column (the sequence) and puts it in the text file roommateMappedReads.txt. Then I would write a python script to count the total number of characters (excluding new line chars), and output that result / 1665.  

## Citations
1. Cermak, Vladimir. “WebDSV.” WebDSV - Free Online DNA Sequence Editor, www.molbiotools.com/WebDSV/index.html.  

2. “Influenza (Flu).” Centers for Disease Control and Prevention, Centers for Disease Control and Prevention, 30 Mar. 2018, www.cdc.gov/flu/.

3. “Influenza A Virus (A/USA/RVD1_H3/2011(H3N2)) Segment 4 Hemagglutinin ( - Nucleotide - NCBI.” National Center for Biotechnology Information, U.S. National Library of Medicine, www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta.  

4. Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research DOI: 10.1101/gr.129684.111  

5. Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]  

6. Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]  

7. Muñoz, Enrique T., and Michael W. Deem. “Epitope Analysis for Influenza Vaccine Design.” Vaccine 23.9 (2005): 1144–1148. PMC. Web. 14 Apr. 2018.  

8. Robasky K, et al. The role of replicates for error mitigation in next-generation sequencing, Nat. Rev. Genet. , 2014, vol. 15 (pg. 56-62)  
