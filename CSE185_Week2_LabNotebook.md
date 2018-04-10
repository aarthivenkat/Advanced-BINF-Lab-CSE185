# Why did I get the flu?
#### Name: Aarthi Venkat
#### Date: Apr 10 2018

Your roommate has a virus with a profile closely matching A/Hong Kong/4801/2014 (H3N2). This flu strain was covered in this season vaccine. Thus, because I got my flu vaccine, I posit that a small portion of the virus population mutated and evolved while replicating inside my roommate. We perform a single-end sequencing run to analyze the sequence of the virus in my roommate.

## 1. Inspecting the data from your roommate  

**How many reads are in roommate.fastq? (iClicker)**  
`wc -l roommate.fastq` gives us the linecount of the file. 

1146956 roommate.fastq  

Because fastq files have 4 lines per read, there are 286739 reads.  

**Look at the first 20 lines and answer (iClicker)**  
`head -n 20 roommate.fastq`  

@SRR1705889.1 1 length=151
TATTAACCATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTGCCAAACGGAGCGATAGTGAAAACAATCACGAATGA
+SRR1705889.1 1 length=151
?????BBBDDDDDDDDGGGGGGIIIHHIFFFHHHHHHIIHIIIFHHGHDFFHGAEHHHIHIIHHHHHIHFHIHIBFHHIHFFHFHHBE>DHHHHDHHHHCDFF;ADEHDFFHHG?BBBFDFFDDEGGB6AC?A>ACGACEE-CEEE8CCEC
@SRR1705889.2 2 length=151
TATTAACCATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTGCCAAACGGAACGATAGTGAAAACAATCACGAATGA
+SRR1705889.2 2 length=151
????ABBBDDDDDDDDGGGGGGIIHHIIIIIIIHHHIIIIIIIIIIIIGHHIHHHHIHHGHFHIIIIIIGHIIIIIIIIIIIIIIHHHHHIHHHHHHHHHHHHHHHHHHHHHHGFGGGFFGGGGGGGGGGGGACCEGGGGGGEGGGCGGG?
@SRR1705889.3 3 length=131
ATCGTTCCGTTTGGCACTGCATGGTGCCCAAGGCACAGCGTTGCCGTGCTGTTGTCATTTCCAGGAAGTTTTTGAGCGAAAACCAGACATAGAATGTAGCTCAAAGCAATGATAGTCTTCATGGTTAATAG
+SRR1705889.3 3 length=131
??????BBDDDDDDDDGGGGGGIIFHIIIIIIIIIIIIIHHHHHIHEEHIIHHIIIIIIIIIIHIFHHHIIIIHHIIHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGGGGGGGGGGGGGGGGDEGGGGG
@SRR1705889.4 4 length=131
ATCGTTCCGTTTGGCACTGCATGGTGCCCAAGGCACAGCGTTGCCGTGCTGTTGTCATTTCCAGGAAGTTTTTGAGCGAAAACCAGACATAGAATGTAGCTCAAAGCAATGATAGTCTTCATGGTTAATAG
+SRR1705889.4 4 length=131
??????BBDDDDEEEEGGGGGGIIFHIIIIIIIGIIIIHEHHHHIHHHHHHHHHHHIIIIIIIIIHHIGHHIIHHIIHHHHHHHHHHHHHHHHHHHHFHHHHFHHHGFGGGGGGGDDEGGGEGGEGGEEGG
@SRR1705889.5 5 length=108
GTGCCCAAGGCACAGCGTTGCCGTGCTGTTGTCATTTCCAGGAAGTTTTTGAGCGAAAACCAGACATAGAATGTAGCTCAAAGCAATGATAGTCTTCATGGTTAATAG
+SRR1705889.5 5 length=108
?????BB?BB9BBBBBC@CA>CCEE>E;EFF7CGDCFAA9CEAFFEEEFHC>EDEHHHC=>+AEFFHHFGGHHDGHHHGHHHH?DDD=DGHFHHF.7D@C..7CD,C,

roommate.fastq contains reads of multiple lengths, so we run the following command to sorted unique read lengths (col2) and the count of each read length (col1).  

`cat roommate.fastq | awk 'NR%4==0 {print length}' | sort -n | uniq -c`  

     74 35  
     16 36  
     24 37  
     36 38  
     30 39  
     37 40  
     ...
   3152 146  
  13708 147  
   5400 148  
  15570 149  
  45169 150  
 187237 151  

So, the maximum read length is 151 bp. **Answer iClicker again.**  

## 2. Alignment of roommate data to reference sequence  

First, we download the reference sequence with NCBI id number KF848938.  

`efetch -db nucleotide -id KF848938.1 -format fasta > KF848938.1.fasta`  

Then, we index the reference file and align the sample viral data to the reference sequence. Because we learned how to do this procedure in the last lab, we can simply pipe the output of bwa mem into samtools and create the bam file.  

`bwa index KF848938.1.fasta`  
`bwa mem KF848938.1.fasta /home/linux/ieng6/cs185s/public/week2/roommate.fastq | samtools view -S -b | samtools sort > roommate.bam`  

The samtools view -f 4 command will allow us to extract all the unmapped reads from the bam file, and we can go on to count these reads and calculate the number mapped.  

`samtools view -f 4 roommate.bam | wc-l`  
3430  


So, if there were initially 286739 reads, and 3430 are unmapped, then 283309 reads were mapped.  

We finally index the bam file using samtools index.  

`samtools index roommate.bam`  


