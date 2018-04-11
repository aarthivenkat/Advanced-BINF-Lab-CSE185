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
For each fastq file f, we `wc -l f` and divide the result by 4.  

SRR1705858.fastq    256586 reads
SRR1705859.fastq    233327 reads
SRR1705860.fastq    249964 reads  

We should calculate a rough estimate of the coverage. Coverage is the total number of base pairs read in the FASTQ files divided by the total number of base pairs in the reference sequence.  

`wc -l KF848938.1.fasta`  
25  
25 - 2 (header and last line) = 23 * 70 + 55 (last line)  

**???** (249964 reads * 151 bp/read) / (23*70 + 55 bp) = 22669x coverage  

#### To align the control fastq files to the reference and create the mpileup file for VarScan,   
```
for x in SRR1705858 SRR1705859 SRR1705860
do
  echo "Aligning $x"
  bwa mem KF848938.1.fasta /home/linux/ieng6/cs185s/public/week2/$x.fastq | samtools view -S -b | samtools sort > $x.bam
  samtools index $x.bam
  samtools mpileup -d 1000000 -f KF848938.1.fasta $x.bam > $x.mpileup`
done
```  

#### Next, we run VarScan with minimum variant frequency 0.1%, only outputting variants and formatting in VCF.  
```
for x in SRR1705858 SRR1705859 SRR1705860
do
  echo "VarScan on $x"
  java -jar /home/linux/ieng6/cs185s/public/tools/VarScan.jar mpileup2snp $x.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > $x.vcf
done
```  

#### To parse the VCF file and get the pos, refbase, altbase, and freq,  
```
for x in SRR1705858 SRR1705859 SRR1705860
do  
  echo "Parsing $x"
  cat $x.vcf | grep -v "^#" | awk '{print $2, $4, $5, $10}' | awk -F '[ :]' '{print $1, $2, $3, $10}'
done
```

#### SRR1705858 VARIANTS  

38 T C 0.66%    
54 T C 0.3%     
72 A G 0.3%     
95 A G 0.24%    
117 C T 0.3%    
165 T C 0.24%   
183 A G 0.3%    
216 A G 0.22%   
218 A G 0.28%   
222 T C 0.26%   
235 T C 0.25%   
254 A G 0.25%   
276 A G 0.22%   
297 T C 0.2%    
328 T C 0.2%    
340 T C 0.23%   
356 A G 0.22%   
370 A G 0.21%   
389 T C 0.23%   
409 T C 0.22%   
414 T C 0.28%   
421 A G 0.18%   
426 A G 0.19%   
463 A G 0.19%   
516 A G 0.2%    
566 A G 0.22%   
595 G T 0.34%   
597 A G 0.17%   
660 A G 0.2%    
670 A G 0.29%   
691 A G 0.23%   
722 A G 0.23%   
744 A G 0.21%   
774 T C 0.3%    
859 A G 0.27%   
915 T C 0.26%   
987 A G 0.22%   
1008 T G 0.27%  
1031 A G 0.28%  
1043 A G 0.24%  
1056 T C 0.2%   
1086 A G 0.33%  
1089 A G 0.22%  
1213 A G 0.24%  
1260 A C 0.3%   
1264 T C 0.26%  
1280 T C 0.25%  
1281 T C 0.22%  
1286 T C 0.2%   
1339 T C 0.41%  
1358 A G 0.26%  
1398 T C 0.2%   
1421 A G 0.31%  
1460 A G 0.34%  
1482 A G 0.24%  
1580 T C 0.25%  
1591 T C 0.29%  

### SRR1705859 VARIANTS  
44 T C 0.47%    
158 A G 0.24%   
165 T C 0.27%   
183 A G 0.22%   
193 A G 0.22%   
216 A G 0.24%   
218 A G 0.29%   
222 T C 0.25%   
254 A G 0.19%   
276 A G 0.24%   
319 T C 0.23%   
340 T C 0.21%   
356 A G 0.24%   
370 A G 0.21%   
398 A G 0.22%   
403 A G 0.19%   
409 T C 0.19%   
414 T C 0.22%   
421 A G 0.18%   
463 A G 0.19%   
499 A G 0.21%   
516 A G 0.2%    
548 A G 0.19%   
591 A G 0.19%   
607 A G 0.18%   
660 A G 0.27%   
670 A G 0.28%   
691 A G 0.23%   
722 A G 0.23%   
744 A G 0.25%   
793 A G 0.17%   
859 A G 0.29%   
898 A G 0.2%    
915 T C 0.21%   
987 A G 0.22%   
1031 A G 0.28%  
1056 T C 0.19%  
1086 A G 0.21%  
1100 T C 0.21%  
1213 A G 0.22%  
1264 T C 0.21%  
1280 T C 0.24%  
1358 A G 0.25%  
1366 A G 0.22%  
1398 T C 0.23%  
1421 A G 0.24%  
1460 A G 0.37%  
1482 A G 0.25%  
1517 A G 0.24%  
1520 T C 0.27%  
1600 T C 0.35%  
1604 T C 0.31%  

### SRR1705860 VARIANTS  
38 T C 0.7%         
44 T C 0.5%         
95 A G 0.24%        
105 A G 0.25%       
133 A G 0.22%       
158 A G 0.26%       
165 T C 0.25%       
183 A G 0.23%       
199 A G 0.19%       
216 A G 0.24%       
218 A G 0.23%       
222 T C 0.3%        
228 T C 0.19%       
230 A G 0.19%       
235 T C 0.25%       
254 A G 0.23%       
271 A G 0.21%       
276 A G 0.33%       
297 T C 0.23%       
319 T C 0.21%       
340 T C 0.21%       
356 A G 0.21%       
370 A G 0.22%       
389 T C 0.2%        
409 T C 0.19%       
414 T C 0.3%        
421 A G 0.21%       
463 A G 0.2%        
499 A G 0.19%       
566 A G 0.24%       
597 A G 0.18%       
607 A G 0.2%        
660 A G 0.28%       
670 A G 0.33%       
691 A G 0.23%       
722 A G 0.25%       
744 A G 0.22%       
759 T C 0.19%       
859 A G 0.25%       
915 T C 0.27%       
987 A G 0.22%       
1031 A G 0.26%      
1043 A G 0.21%      
1056 T C 0.2%       
1086 A G 0.3%       
1089 A G 0.22%      
1105 A G 0.22%      
1209 A G 0.27%      
1213 A G 0.24%      
1264 T C 0.27%      
1280 T C 0.25%      
1281 T C 0.21%      
1301 A G 0.22%      
1358 A G 0.29%      
1366 A G 0.21%      
1398 T C 0.23%      
1421 A G 0.37%      
1460 A G 0.26%      
1482 A G 0.23%      
1580 T C 0.27%      
1604 T C 0.3%      

## 5. Compare the control results to your roommate's results    
Control|AvgFreq|StdevFreq  
--|--|--  
SRR1705858|0.256%|0.0717%  
SRR1705859|0.237%|0.0524%  
SRR1705860|0.250%|0.0780%  

**iClicker**  

Did VarScan report rare mutations in your roommate's file with frequencies that are more than 3 standard deviations from the averages in the reference files?  

Yes. Using the excel spreadsheet with the controls and roommate variants (excluding common variants with frequency >95%), I wrote the nested if statement  

`=IF(FREQ>(0.256%+(3*0.0717%)),">3std",IF(FREQ>(0.237%+(3*0.0524%)),">3std",IF(FREQ>(0.25%+(3*0.078%)),">3std","<3std")))`  

**Question: should we check if 3std > than ALL averages, or AT LEAST ONE of averages? Or do we average the averages/std?**  

Where FREQ refers to the variant allele frequency for each rare variant. VarScan reported the following variants as >3std from any one of the averages of the reference files.  

Position|RefAllele|AltAllele|VariantAlleleFreq  
---|---|---|---  
38|T|C|0.45%  
495|C|T|1.04%  
910|G|A|0.73%  
1293|G|A|61.82%
1521|G|A|1.12%  

We'll turn again to WebDSV to understand how the amino acids change for these variants. We want to record the original codon, mutated codon, original amino acid, its position in the protein, and the mutated amino acid. Then, record whether the change is synonymous or missense. Example: A72G ACA>ACG Thr24Thr synonymous.  

Variant|Codon Change|AA Change/Position|Mutation Type  
---|---|---|---  
T38C|CTG>CCG|Leu13Pro|missense  
C495T|AAC>AAT|Asn165Asn|synonymous  
G910A|GCC>ACC|Ala304Thr|missense
G1293A|CTG>CTA|Leu431Leu|synonymous
G1521A|CTG>CTA|Leu507Leu|synonymous  

**Discussion question: Are there any positions reported by VarScan in all 3 of the reference sequences? You could, in principle, also calculate the average and standard deviation between the 3 reference replicates for one position at a time. Which kind of average and standard deviation do you think is better for error correction?**  

The following positions are reported by VarScan for all 3 of the reference sequences. Because this list contains positions that consistently poorly sequenced, it may be wiser for error correction to compute the average and standard deviation of variant allele frequency at only these positions.   
165
183
216
218
222
254
276
340
356
370
409
414
421
463
660
670
691
722
744
859
915
987
1031
1056
1086
1213
1264
1280
1358
1398
1421
1460
1482  

## 6. Epitope Mapping

**Use the epitope locations listed in Munoz et al (listed under reading for the lab) to determine if any of the high confidence (> 3 std deviations away from reference error rate) mutations from your roommate’s flu infection are located in an epitope region of hemeagglutinin. (Epitopes are the parts of the protein structure recognized by antibodies). If so, list which epitope regions are mutated**  

**QUESTION. Is this correct?**  

The residues from my roommate's flu infection are 13, 165, 304, 431, and 507.  
Residues 13, 431 and 507 are **not** located in an epitope region.  
Reside 165 is located in **Epitope B** of hemeagglutinin.  
Reside 304 is located in **Epitope C** of hemeaglutinin.  


## 7. References 
[Flu Strains 2017](https://www.cdc.gov/flu/about/season/flu-season-2017-2018.htm)  
[KF848938.1 FASTA](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta)  
[BWA](http://bio-bwa.sourceforge.net/)  
[SamTools](http://samtools.sourceforge.net/)  
[VarScan](http://varscan.sourceforge.net/)  
[WebDSV](http://www.molbiotools.com/WebDSV/)
