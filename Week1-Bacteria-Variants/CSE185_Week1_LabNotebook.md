# CSE 185 Lab Notebook - Week 1

#### Name: Aarthi Venkat
#### Date: April 3, 2018

## 1. UNIX Navigation: Determining File Sizes

`ls -lh`

-rwxr-xr-x 1 cs185s ieng6_cs185s 4.5M Mar 12 13:46 NC_000913.3.fasta  
-rwxr-xr-x 1 cs185s ieng6_cs185s 1.2G Mar 12 13:46 amp_res_1.fastq  
-rwxr-xr-x 1 cs185s ieng6_cs185s 1.2G Mar 12 13:46 amp_res_2.fastq  

The file sizes are **4.5M**, **1.2G**, **1.2G**, for the fasta file and the fastq files, respectively.  

## 2. Raw Data Inspection: Determining Quality Manually

`head -n 20 amp_res_1.fastq`  

Each read has 4 lines of information. The first line starts with @ and contains identifiers and information about the read. The next line is the actual sequence. The third line starts with a + and repeats identifier information. The fourth line contains the quality string, whree ASCII characters encode the quality score for each base. Below is the 3rd read of amp_res_1.fastq.  

@SRR1363232.3 GWZHISEQ01:153:C1W31ACXX:5:1101:3103:2142 length=101
CTAAGATACACTCAACTGATATAGCNTCTCTCTACTACTATNGNNNNNTAGTCCACTTTCATAATGAAAAAATTGAAGAGGCAAGGATTTGTATAGACAAA
+SRR1363232.3 GWZHISEQ01:153:C1W31ACXX:5:1101:3103:2142 length=101
@@@ADFDFDBDDFHGEGIEIG?C@E#2AECGE@EGIIIJEI#0#####0077BF=C@BFC@GGHGGCHGGHDFDFFEEA;A@@??=<>@A@DDEC@CDCC9  

The quality symbol for the read's first base is @, which corresponds to a **quality score of 31** according to [the given resources](http://drive5.com/usearch/manual/quality_score.html).  

## 3. Raw Data Inspection: Calculating Number of Reads  

`wc -l amp_res_1.fastq`  

There are 14214080 lines in the document, and because we know each read as 4 lines, **there are 3553520 reads in the file**.  
Similarly, the amp_res_2.fastq file has 3553520 reads. The fasta file is the line count - 1 (to account for the header), which results in 66310 reads.

## 4. Raw Data Inspection: Coverage  

Given that the length of the genome for this strain is 4641652 base pairs, we determine the average coverage by multiplying the number of reads of each fastq file by the length of each read. Knowing there are 3553520 reads in both files, and each read is 101 bp, the total base pairs we have read is 717811040, and **the coverage with respect to genome is 155x**.  

## 5. Running FastQC: High-Throughput Sequence QC Analysis Tool  

FastQC reads a set of sequence files and produces from each one a quality control report consisting of a number of different modules, each one of which will help to identify a different potential type of problem in your data.  

#### Running FastQC given two FASTQ files  
`fastqc -o . /home/linux/ieng6/cs185s/public/week1/amp_res_1.fastq /home/linux/ieng6/cs185s/public/week1/amp_res_2.fastq`  

#### Copying HTML output to Desktop  
After maneuvering to local Desktop:  
`scp avenkat@ieng6.ucsd.edu:/home/linux/ieng6/cs185s/avenkat/week1/amp_res_*_fastqc.html .`  

#### Do the basic statistics match what you calculated for the number of reads last time?  
Yes, the basic statistics match.

## 6. Analyzing FastQC failures  

#### Red Circles
For amp_res_1.fastq, there is a red circle (indicating irregular results) for Per Base Sequence Quality and Per Tile Sequence Quality.  
For amp_res_2.fastq, there is a red circle for Per Base Sequence Quality and Per Base Sequence Content.  

#### Per Base Sequence Quality  
According to ![FastQC Documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/), Per Base Sequence Quality refers to the analysis of the quality per base as a box plot.  

*This will result in failure when, for any base, the lower quartile of the box plot is <5 or the median is <20.*

This error can occur for one of two reasons:  
1. General degradation over duration of long runs  
2. Short loss of quality early in the run  

The remedy for the first issue is *quality trimming*, which truncates reads based on average quality. The remedy for the second issue is masking bases during subsequent assembly.  

#### Per Tile Sequence Quality  
This quality measure allows you to determine if loss in quality is associated with only one part of the flowcell.  The module will output a failure if any tile indicates a mean Phred score >5 and <mean for that base across all tiles. This error generally occurs because the flowcell was overloaded, and as such events are not confined to a certain area of the flowcell.  

#### Per Base Sequence Content  
This measure indicates the proportion of each base at each position in a file. A failure occures if the difference between any of the four bases is >20% at any position. This could have occured for a number of reasons relating to overrepresentation of certain sequences and biases in preparation.  

#### What should we do about things fastqc identified as unusual (iClicker)?  
We must analyze why these unusual outputs occur before continuing with subsequent analyses. If there are significant issues with the data at hand, we cannot produce valuable output downstream.  

## 7. Filtering the reads with Quality Trimming  
The trimming program we will use is Sickle. According to the lab instructions, Sickle uses a small window to slide across the sequence, trimming until the average quality is above a specified threshold, and cutting off when average quality drops below.  

`sickle pe -f /home/linux/ieng6/cs185s/public/week1/amp_res_1.fastq -r /home/linux/ieng6/cs185s/public/week1/amp_res_2.fastq -t sanger -o trimpair1.fastq -p trimpair2.fastq -s singletons.fastq`  

PE forward file: /home/linux/ieng6/cs185s/public/week1/amp_res_1.fastq  
PE reverse file: /home/linux/ieng6/cs185s/public/week1/amp_res_2.fastq  

Total input FastQ records: 7107040 (3553520 pairs)  

FastQ paired records kept: 6904494 (3452247 pairs)  
FastQ single records kept: 98694 (from PE1: 94870, from PE2: 3824)  
FastQ paired records discarded: 5158 (2579 pairs)  
FastQ single records discarded: 98694 (from PE1: 3824, from PE2: 94870)  

#### Manually Checking Number of Trimmed Reads  
`wc -l trimpair1.fastq`  
13808988 trimpair1.fastq  
This line count divided by 4 results in 3452247 reads.  

`wc -l trimpair2.fastq`  
13808988 trimpair1.fastq  
This line count divided by 4 results in 3452247 reads.  

This matches the above output for number of total reads and number of pairs.  

#### Rerunning FastQC on Trimmed Reads  
`fastqc -o . trimpair1.fastq trimpair2.fastq`

#### Copying HTML output to Desktop  
After maneuvering to Desktop  
`scp avenkat@ieng6.ucsd.edu:/home/linux/ieng6/cs185s/avenkat/week1/trimpair*_fastqc.html .`  

Now, the quality per base does not go below a quality threshold of 20.  

#### Threshold Incresed to 30  
`sickle pe -q 30 -f /home/linux/ieng6/cs185s/public/week1/amp_res_1.fastq -r /home/linux/ieng6/cs185s/public/week1/amp_res_2.fastq -t sanger -o trimpair1_30.fastq -p trimpair2_30.fastq -s singletons_30.fastq`

PE forward file: /home/linux/ieng6/cs185s/public/week1/amp_res_1.fastq
PE reverse file: /home/linux/ieng6/cs185s/public/week1/amp_res_2.fastq

Total input FastQ records: 7107040 (3553520 pairs)

FastQ paired records kept: 6012550 (3006275 pairs)
FastQ single records kept: 375230 (from PE1: 283732, from PE2: 91498)
FastQ paired records discarded: 344030 (172015 pairs)
FastQ single records discarded: 375230 (from PE1: 91498, from PE2: 283732)

#### How many paired reads did sickle keep after using a quality score threshold of 30 for trimming (iClicker)?  
As determined from the above output, Sickle kept **3006275 pairs** with the higher threshold.

## 8. Read Alignment with BWA-MEM  

BWA-MEM uses a Burrows-Wheeler transform for data compression to find maximum exact match within a read to the reference.  

First, we copy the reference genome into a local folder.  
`cp /home/linux/ieng6/cs185s/public/week1/NC_000913.3.fasta .`

Then, we index sequences in the FASTA format.  
`bwa index NC_000913.3.fasta`  

The following output files are created:  
NC_000913.3.fasta.amb  
NC_000913.3.fasta.ann  
NC_000913.3.fasta.bwt  
NC_000913.3.fasta.pac  
NC_000913.3.fasta.sa  

To run the algorithm,  
`bwa mem NC_000913.3.fasta trimpair1.fastq trimpair2.fastq > trimpaired.sam`  

## 9. Samtools and SAM  

`head -n 5 trimpaired.sam | less -S`  
**How many separate lines are in the output from head? How many reads (iClicker)?**  
There are 5 lines in the output and 3 reads.  

Using samtools to read and parse sam files  

`samtools flagstat trimpaired.sam`  

6906203 + 0 in total (QC-passed reads + QC-failed reads)  
0 + 0 secondary  
1709 + 0 supplementary  
0 + 0 duplicates  
6899109 + 0 mapped (99.90% : N/A)  
6904494 + 0 paired in sequencing  
3452247 + 0 read1  
3452247 + 0 read2  
6881172 + 0 properly paired (99.66% : N/A)  
6891682 + 0 with itself and mate mapped  
5718 + 0 singletons (0.08% : N/A)  
0 + 0 with mate mapped to a different chr  
0 + 0 with mate mapped to a different chr (mapQ>=5)  

The count paired in sequencing is 6904494, which does not match the in total of 6906203. 99.90%, or 6899109 reads, are mapped.

`samtools view -S -b trimpaired.sam > trimpaired.bam`
`samtools sort trimpaired.bam > sortedtrimpaired.bam`
`samtools index sortedtrimpaired.bam`

To visualize alignment of the reads, we can run  
`samtools tview sortedtrimpaired.bam NC_000913.3.fasta`  

**How do you interpret the A at position 46 in the alignment (iClicker)?**  
This read contains a sequencing error. Because this A only exists in one read, it is likely that that read contained a sequencing error rather than a mutation, as if there were a mutation, it would exist in multiple reads.  

## 10. Make a pileup, call variants  

We would like go through the BAM file and, for each position in the reference genome, see how many reads have a mutation at the same position. We can't use tview or the SAM file directly for this purpose, so we resort to mpileup through samtools.  

`samtools mpileup -f NC_000913.3.fasta sortedtrimpaired.bam > my.mpileup`

Using `head -n 100`, we can see that this is much easier to read. However, to extract interesting information we will downloard VarScan.  

`curl -L https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download > VarScan.jar`  

Now we need to decide the minimum variant frequency to accurately identify true mutations in our antibiotic strain. **Where do we set the minimum variant frequency (iClicker)?**  
We can set the minimum threshold to 70%. As the clicker answer states, there shouldn't be more than one kind of bacteria in my sample, but I'm wiling to tolerate some sequencing errors.  

To run this program with a threshold of 70%,  
`java -jar VarScan.jar mpileup2snp my.mpileup --min-var-freq 0.7 --variants --output-vcf 1 > VarScan.vcf`  

Only SNPs will be reported  
Warning: No p-value threshold provided, so p-values will not be calculated  
Min coverage:	8  
Min reads2:	2  
Min var freq:	0.7  
Min avg qual:	15  
P-value thresh:	0.01  
Reading input from my.mpileup  
4641643 bases in pileup file  
9 variant positions (6 SNP, 3 indel)  
0 were failed by the strand-filter  
6 variant positions reported (6 SNP, 0 indel)  

VarScan.vcf has the following information indicating mutations in the E.coli strain. I have recorded the position, reference, alternative base, and variant allele frequency in each column in the table below.

POS | REF BASE | ALT BASE | FREQ
---|---|---|---
92439|G|A|99.32%
803662|C|A|100%
852762|A|G|100%
1905761|A|G|100%
3535147|G|A|100%
4390754|G|T|100%

## 11. Variant effect prediction  

In order to use the webtool Ensemble Variant Effect Predictor, we change the VarScan.vcf to account for the chromosome format in the first field.  

`awk '{if (NR>24) $1="Chromosome"; print}' VarScan.vcf > mymodVarScan.vcf`  

Moving this VCF file to my local Desktop  

`scp avenkat@ieng6.ucsd.edu:/home/linux/ieng6/cs185s/avenkat/week1/mymodVarScan.vcf .`  

Using the Variant Effect Predictor tool for E.coli, and returning results for variants in coding regions only, we get the following results:  


POS | OCCURS IN GENE? | TYPE | GENE NAME | AA SUBSTITUTION
---|---|---|---|---
92439|yes|missense|b0084|D to N
803662|yes|missense|b0771|L to I
852762|no|N/A|N/A|N/A
1905761|yes|missense|b1821|G to D
3535147|yes|missense|b3404|V to G
4390754|yes|synonymous|b4161|A to A  
