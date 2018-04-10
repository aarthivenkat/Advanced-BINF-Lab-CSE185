# Why did I get the flu?
#### Name: Aarthi Venkat
#### Date: Apr 10 2018

Your roommate has a virus with a profile closely matching A/Hong Kong/4801/2014 (H3N2). This flu strain was covered in this season vaccine. Thus, because I got my flu vaccine, I posit that a small portion of the virus population mutated and evolved while replicating inside my roommate. We perform a single-end sequencing run to analyze the sequence of the virus in my roommate.

## 1. Inspecting the data from your roommate  

**How many reads are in roommate.fastq? (iClicker)**  
`wc -l roommate.fastq` gives us the linecount of the file. 

1146956 roommate.fastq  

Because fastq files have 4 lines per read, there are 286739 reads.  

**Based on the first 20 lines of the fastq file, how many cycles were probably carried out during the sequencing run? (iClicker)**  
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

roommate.fastq contains reads of multiple lengths, and the reads are not sorted by length in any definitive way, so it is currently impossible to tell for certain how many cycles were carried out. To determine the number of cycles, we run the following command to sorted unique read lengths (col2) and the count of each read length (col1).  

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

So, the maximum read length is 151 bp, so unless every read was trimmed, it is likely that the number of cycles carried out during the sequencing run is 151.  

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

## 3. Look for common variants with VarScan  

We know want to make an mpileup of the bam alignment file using samtools. The default depth limit is 8000 calls, but because are looking into rare variants, we are increasing the depth limit to 1000000.  

`samtools mpileup -d 1000000 -f KF848938.1.fasta roommate.bam > roommate.mpileup`  

Now we can run VarScan on mpileup. We initially set a high minimum variant frequency cut-off so we look only at common mutants (those present in 95% of viral DNA molecules).  

`java -jar /home/linux/ieng6/cs185s/public/tools/VarScan.jar mpileup2snp roommate.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > roommate.vcf`  


We are interested in files 2,4 and 5 in the VCF file (excluding the header), so we output only these fields using the following script:  

`cat roommate.vcf | grep -v "^#" | awk '{print $2, $4, $5}'`  

It is clear now that the following 7 variants reported back.  

Reference Base | Position | Mutant Base  
---|---|---
A|72|G
C|117|T
G|595|T
T|774|C
T|1008|G
A|1260|C
T|1339|C


We are interested in finding how the codons and amino acids changed with the mutation, so we copy the reference FASTA sequence into online sequence editor WebDSV.  

For each variant, we want to record the original codon, mutated codon, original AA, position in protein, mutated AA, and whether the change is synonymous or missense.  

Variant | Original Codon > Mutated Codon | Original AA + Position in Protein + Mutated AA | Type of Change
---|---|---|---
A72G | ACA > ACG| Thr24Thr | synonymous
C117T | GCC > GCT | Ala39Ala | synonymous
G595T | GCA > TCA | Ala199Ser | missense
T774C | TTT > TTC | Phe258Phe | synonymous
T1008G | GCT > GCG | Ala336Ala | synonymous
A1260C | CTA > CTC | Leu420Leu | synonymous
T1339C | TTG > CTG | Leu447Leu | synonymous

**How many common missense viariants did you find in the sample from your roommate? (iClicker)** 
Of the six variants, I found one missense variant.  

## 4. Look for rare variants with VarScan  

Now, we repeat our analysis of the mpileup file, but we set the minimum var frequency to 0.1% to search for rare variants.  

`java -jar /home/linux/ieng6/cs185s/public/tools/VarScan.jar mpileup2snp roommate.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > roommate_rare.vcf`  

We want to look at the VCF file to understand the frequency of each of the variants.  
`cat roommate_rare.vcf | grep -v "^#" | awk '{print $2, $4, $5, $10}'`  

38 T C 0/1:22:3113:3111:3097:14:0.45%:6.3047E-3:35:31:2268:829:12:2  
72 A G 1/1:255:9835:9820:2:9817:99.97%:0E0:31:36:2:0:6811:3006  
117 C T 1/1:255:13852:13732:12:13719:99.91%:0E0:33:37:8:4:9670:4049  
...

Evidently, the frequency 8 is buried in the last field. So, we pipe this output into more awk command to process the output futher.  

`cat roommate_rare.vcf | grep -v "^#" | awk '{print $2, $4, $5, $10}' | awk -F '[ :]' '{print $1, $2, $3, $10}'`  

The results of this command are as follows:  

38 T C 0.45%  
72 A G 99.97%  
117 C T 99.91%  
216 A G 0.18%  
218 A G 0.2%  
254 A G 0.28%  
276 A G 0.37%  
295 C T 0.23%  
319 T C 0.22%  
409 T C 0.26%  
495 C T 1.04%  
524 A G 0.19%  
595 G T 99.94%  
691 A G 0.17%  
722 A G 0.19%  
774 T C 99.96%  
910 G A 0.73%  
915 T C 0.27%  
987 A G 0.25%  
1008 T G 99.9%  
1043 A T 0.19%  
1086 A G 0.26%  
1100 T C 0.2%  
1260 A C 99.9%  
1293 G A 61.82%  
1339 T C 99.97%  
1460 A G 0.23%  
1473 C T 0.23%  
1517 A G 0.22%  
1521 G A 1.12%  
1604 T C 0.25%  

In order to determine which are rare variants and which are sequencing errors, we perform control experiments with three samples.  

SRR1705858.fastq    256586 reads
SRR1705859.fastq    233327 reads
SRR1705860.fastq    249964 reads

## 5. References 
https://www.cdc.gov/flu/about/season/flu-season-2017-2018.htm  
http://www.molbiotools.com/WebDSV/index.html  
https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta  
[BWA](http://bio-bwa.sourceforge.net/)  
[SamTools](http://samtools.sourceforge.net/)  
