# CSE 185 Lab Notebook - Week 3  

## Inspect the data and run fastqc  

#### Number of reads in each raw data file  
`wc -l`/ 4  

file|number of reads  
---|---  
frag_1.fastq|1557594  
frag_2.fastq|1557594  
short_jump_1.fastq|1111432  
short_jump_2.fastq|1111432  

#### Running fastqc on all four files
`fastqc frag_1.fastq frag_2.fastq`  
`fastqc short_jump_1.fastq short_jump_2.fastq`  
HTML output moved to labreport/  
fastq and zip files ignored in .gitignore file  

## Use jellyfish to make k-mer histograms for the "frag" library  

Jellyfish is a kmer counting program that counts the frequency of all possible kmers of a given length in the data. With the below command, we specify the length 31, ignore directionality, estimate a large size of the hashtable, and set the name of the file to be the kmer length.  

`jellyfish count -m 31 -s 10000000 -o 31 -C frag_1.fastq`  

To make the histogram file,  
`jellyfish histo 31 > 31.histo`  
column one: list of bins (number of times a kmer occurs)  
column two: count for the number of kmers in that category  

## Plot the kmer histograms in python  

```
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.figure()
r31 = pd.read_csv("31.histo", sep=" ", names=["kmercount", "number"])
plt.bar(r31.iloc[1:100,]["kmercount"], r31.iloc[1:100,]["number"]);
plt.savefig("frag_1_31.pdf")
```
**iClicker #1**  
After moving the histogram to github, I recognize that the very strong peak in the beginning of the plot is likely sequencing errors.  

**iClicker #2**  
The bin that corresponds to the *initial valley point* is bin 5, with 34039 kmers fitting that category.  

## K-mer based error correction  
`vim filelist`  
Add frag_1.fastq and frag_2.fastq to the filelist  

`KmerFreq_HA -k 27 -L 101 -i 10000000 -p corrected -l filelist` 

According to the tutorial, this command is generating the frequency of every k-mer in our data, then creating a hash-table to make them easier to access.  


`Corrector_HA -k 27 -Q 33 -o 3 -l 5 corrected.freq.gz filelist`  
This command corrects the data based on the KmerFreq analysis.  

```
jellyfish count -m 31 -s 10000000 -o 31corrected -C frag_1.fastq.cor.pair_1.fq  
jellyfish histo 31corrected > 31corrected.histo  
```
We will see how the distribution of kmers is different from the corrected data by running jellyfish once again.  

It looks better!

## Collect data on corrected reads, calculate genome size  

#### Number of reads in each corrected data file
`wc -l`/ 4
file|number of reads  
---|---  
frag_1.fastq.cor.pair_1.fq|1054393  
frag_2.fastq.cor.pair_2.fq|1054393  


#### New valley point and peak  
The new valley point is the bin 2, with frequency 15515. 
The peak of the graph is the bin 14, with frequency 142225.  

#### Average read length of corrected data files  
`awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' x.fq` for x in file column below.  

file|average read length  
---|---  
frag_1.fastq.cor.pair_1.fq|95.631  
frag_2.fastq.cor.pair_2.fq|93.3334  


#### Total number of bases of corrected data files  
`awk 'NR%4==2{sum+=length($0)}END{print sum}' frag_1.fastq.cor.pair_1.fq`  

file|total number of base pairs  
---|---  
frag_1.fastq.cor.pair_1.fq|100832704  
frag_2.fastq.cor.pair_2.fq|98410110  

Using the formulas to calculate genome size N = (M*L)/(L-K+1) Genome_size = T/N  

file|genome size  
---|---
frag_1.fastq.cor.pair_1.fq|5244176.658  

#### Optional extra credit  
The actual length of the genome is around 3 million bp - I got around 5 million. I may have underestimated the amount of coverage by undercounting the amount of overlap between reads.    

## Assemble reads with minia  

Assembling reads for kmer size 31,  
`minia -in correctedlist -kmer-size 31 -abundance-min 2 -out 31minia_out`  

`wc -l 31minia_out.contigs.fa` outputs there are 1398 contigs, so `699` contigs (2 lines per sequence).  

`cat minia_out.contigs.fa | awk 'NR%2==0{print length}' | sort -n | head` sorts contig lengths from smallest to largest, so the shortest contig is the first value of this output, which is `63`.  

`cat minia_out.contigs.fa | awk 'NR%2==0{print length}' | sort -rn | head` sorts contig lengths from largest to smallest, so the longest contig is the first value of this output, which is `80334`.  

Using QUAST with the contig file, the N50 value results in `25284`.  

## Use short jump library to scaffold minia contigs  

In sspace_library, we specify  
`Lib1 short_jump_1.fastq short_jump_2.fastq 3500 0.5 RF`  

where 3500 is the insertsize, 0.5 is insertsizeerror, and RF is orientation of an outie.  

Run SSPACE,  
`SSPACE_Basic_v2.0.pl -l sspace_library -s 31minia_out.contigs.fa -z 100 -v 1 -p 1 -b output_scaffold_31`  

In the summary file, we see the following information after scaffolding on raw short jump files.  
field|value   
---|---  
number of scaffolds|284  
max scaffold length|258872  
N50|128380  

## Clean up original shortjump read files, see if scaffolding improves  
We will now clean up short jump reads with quality trimming and redo scaffolding.  

`sickle pe -f short_jump_1.fastq -r short_jump_2.fastq -o trimmedfile1.fastq -p trimmedfile2.fastq -s trimmedsingles.fastq -t sanger`  

rewriting the library file to have the trimmed files as input, and running SSPACE on the trimmed files, we get the following in the summary file  
field|value
---|---
number of scaffolds|414
max scaffold length|120417
N50|41285  


`sickle pe -f short_jump_2.fastq -r short_jump_1.fastq -o trimmedfile1.fastq -p trimmedfile2.fastq -s trimmedsingles.fastq -t sanger`

Lib1 trimmedfile1.fastq trimmedfile2.fastq 3500 0.5 FR

SSPACE_Basic_v2.0.pl -l sspace_library -s 31minia_out.contigs.fa -z 100 -v 1 -p 1 -b output_scaffold_31_trim  

field|value  
---|---  
number of scaffolds|473  
max scaffold length|80334  
N50|25284  

**Which scaffolding was better? The one with the raw data or trimmed?**  
The one with the raw data has a higher N50 score, and thus is a better scaffolding, although I recognize the N50 score does not always account for quality so much as length.  

## Try to close the gaps in the scaffolds  

`GapCloser -a ./Thursday/output_scaffold_31.final.scaffolds.fasta -o output.fa -l 101 -b gap.config`  

We complete the lab evaluating our assembly with QUAST.  
