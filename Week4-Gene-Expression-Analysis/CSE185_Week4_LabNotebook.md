# CSE 185 Lab Notebook - Week 4
#### 04/24/2018 Aarthi Venkat  

## Data inspection and quality control  

**What read length and how many reads were used for each experiment?**  
Assuming same read length for each read,  
```
x=4; for i in *.fq.gz; do y="$(zcat $i | wc -l)"; z="$(zcat $i | head -n 2 | tail -n 1 | wc -c)"; echo $i $((y / x)) $z; done
```  
File NumReads ReadLength  
FL_Rep1_chr5_1.fq.gz 3179399 51  
FL_Rep1_chr5_2.fq.gz 3179399 51  
FL_Rep2_chr5_1.fq.gz 2965029 51   
FL_Rep2_chr5_2.fq.gz 2965029 51   
HL_Rep1_chr5_1.fq.gz 3932838 51   
HL_Rep1_chr5_2.fq.gz 3932838 51   
HL_Rep2_chr5_1.fq.gz 2811913 51   
HL_Rep2_chr5_2.fq.gz 2811913 51   
MB_Rep1_chr5_1.fq.gz 3254975 51   
MB_Rep1_chr5_2.fq.gz 3254975 51   
MB_Rep2_chr5_1.fq.gz 3413939 51   
MB_Rep2_chr5_2.fq.gz 3413939 51   

**FastQC** 
```
for i in /home/linux/ieng6/cs185s/public/week4/*.fq.gz; do fastqc -o ./labreport/ $i; done
```

## RNA-seq sequence alignment  

**SAM Header documenting command to generate BAM file**   
```
for i in *_chr5.bam; do x="$(samtools view $i -H | grep "^@PG")"; echo $x; done
```  

@PG ID:STAR PN:STAR VN:STAR_2.5.3a CL:STAR --runThreadN 4 --genomeDir mm10STAR --readFilesIn FL_Rep1_1.fastq.gz FL_Rep1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix FL_Rep1 --outSAMtype BAM SortedByCoordinate  
@PG ID:STAR PN:STAR VN:STAR_2.5.3a CL:STAR --runThreadN 4 --genomeDir mm10STAR --readFilesIn HL_Rep1_1.fastq.gz HL_Rep1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix HL_Rep1 --outSAMtype BAM SortedByCoordinate  
@PG ID:STAR PN:STAR VN:STAR_2.5.3a CL:STAR --runThreadN 4 --genomeDir mm10STAR --readFilesIn FL_Rep2_1.fastq.gz FL_Rep2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix FL_Rep2 --outSAMtype BAM SortedByCoordinate  
@PG ID:STAR PN:STAR VN:STAR_2.5.3a CL:STAR --runThreadN 4 --genomeDir mm10STAR --readFilesIn HL_Rep2_1.fastq.gz HL_Rep2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix HL_Rep2 --outSAMtype BAM SortedByCoordinate  
@PG ID:STAR PN:STAR VN:STAR_2.5.3a CL:STAR --runThreadN 4 --genomeDir mm10STAR --readFilesIn MB_Rep1_1.fastq.gz MB_Rep1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix MB_Rep1 --outSAMtype BAM SortedByCoordinate  
@PG ID:STAR PN:STAR VN:STAR_2.5.3a CL:STAR --runThreadN 4 --genomeDir mm10STAR --readFilesIn MB_Rep2_1.fastq.gz MB_Rep2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix MB_Rep2 --outSAMtype BAM SortedByCoordinate  

#### samtools view *_chr5.bam shows "N" in the CIGAR scores. What biological feature do you think the "N"s stand for?  

According to the [SAM Specification](https://samtools.github.io/hts-specs/SAMv1.pdf), the "N" describes a skipped region from the reference, which can signify large skips in introns and signify regions that are not transcribed in expression data.  

## Kallisto  

chmod +x run_kallisto.sh  
./run_kallisto.sh   

## Visualizing data using a genome-browser  

genome build: mm10 (from samtools header @PG above)  
Exploring IGV  
![Mecom](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/Mecom.PNG)  

## Comparing overall expression patterns across datasets  

```
paste FL_Rep1/abundance.tsv FL_Rep2/abundance.tsv | cut -f 5,10 | grep -v tpm | awk '(!($1==0 && $2==0))' | datamash ppearson 1:2
```
**0.99712632784422**  

```
for i in FL_Rep1 FL_Rep2 HL_Rep1 HL_Rep2 MB_Rep1 MB_Rep2; do for j in FL_Rep1 FL_Rep2 HL_Rep1 HL_Rep2 MB_Rep1 MB_Rep2; do k="$(paste $i/abundance.tsv $j/abundance.tsv | cut -f 5,10 | grep -v tpm | awk '(!($1==0 && $2==0))' | datamash ppearson 1:2)"; echo $i $j $k; done; done
```

![correlation analysis](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/correlation.PNG)  

**Which tissues were most similar? Most different? How concordant were the replicates? Are replicates more concordant with each other than with other tissues?**  

FL and HL had a correlation rate of 0.959, whereas FL and MB 0.954 and HL and MB 0.931. Therefore, FL and HL were the most similar, and HL and MB the most different. The replicates were highly concordant, with correlation rate of 0.995. All the replicates were more concordant with each other than with the other tissues.  

## Differential expression analysis  

In an R environment, we run the following given code  
```  
library("sleuth")

# Set up the paths to our kallisto results
sample_id = c("FL_Rep1","FL_Rep2","HL_Rep1","HL_Rep2","MB_Rep1","MB_Rep2")
kal_dirs = file.path(sample_id)

# Load metadata
s2c = read.table(file.path("/home/linux/ieng6/cs185s/public/week4/exp_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c = dplyr::mutate(s2c, path = kal_dirs)

# Create sleuth object
so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# Fit each model and test
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

# Get output, write results to file
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
# Note, you may need to edit the output path below if your $HOME
# directory is not the same as your CSE185 course directory
write.table(sleuth_significant, "/home/linux/ieng6/cs185s/avenkat/week4/sleuth_results.tab", sep="\t", quote=FALSE)
```
**How many significant transcripts are there?**  
`wc -l sleuth_results.tab` -1 for header = 1940 significant transcripts  

**What are the gene names for the top 10 genes?**  
rank|target_id|gene name  
---|---|---  
1|ENSMUST00000075453.8|Rpl21  
2|ENSMUST00000002708.4|Shh  
3|ENSMUST00000031131.10|Uchl1  
4|ENSMUST00000031249.7|Sparcl1  
5|ENSMUST00000040576.9|Parm1  
6|ENSMUST00000056355.8|Nat8l  
7|ENSMUST00000058045.4|Garem2  
8|ENSMUST00000079324.13|Ubl3  
9|ENSMUST00000102553.10|Hmgn2  
10|ENSMUST00000112707.2|Lrrc8b  

**For several top hits, find the gene name, navigate to that gene in IGV, and take screenshots to include in your lab report.**  

### Rp121  
![Rpl21](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Rpl21.PNG)  
### Shh  
![Shh](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Shh.PNG)  
### Uchl1  
![Uchl1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Uchl1.PNG)  
### Sparcl1  
![Sparcl1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Sparcl1.PNG)  
### Parm1  
![Parm1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/Parm1.PNG)  

**First take a look at where these marks (H3K27ac and H3K4me1) are falling. Are they near gene regions? Beginning or ends of genes? Other places? Discuss this in the results section of your lab report.**  

H3K27ac peaks at the promoter of Rbm33, quickly falls off, and does not peak in Shh or other near gene regions. H3K4me1 peaks at the promoter of each gene in the region of interest, and gradually falls off, showing higher expression in intronic regions.  

## Extra Credit 1  
Visualize expression of differntially expressed genes as a heatmap (e.g. with transcripts as rows and samples as columns). Cluster the rows by row and column. Do replicates cluster together? Are there clear clusters of up vs. down regulated genes in each tissue?  

```
paste FL_Rep1/abundance.tsv FL_Rep2/abundance.tsv HL_Rep1/abundance.tsv HL_Rep2/abundance.tsv MB_Rep1/abundance.tsv MB_Rep2/abundance.tsv | cut -f 5,10,15,20,25,30 | grep -v tpm | awk '(!($1==0 && $2==0 && $3==0 && $4==0 && $5==0 && $6==0))' > tpm_matrix.txt
```
Then, we run a python script to create a clustermap, which clusters the matrix by row and column.  
```
import matplotlib; matplotlib.use('Agg')
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

f = open("tpm_matrix.txt", 'r')
m = []

for i,line in enumerate(f):
    col = [float(x) for x in line.strip().split('\t')]
    newcol = []
    for x in col:
        if x > 0: newcol.append(math.log10(x))
        else: newcol.append(x)
    m.append(newcol)

map = np.array(m)
plt.figure(figsize=(50,50))
g = sns.clustermap(map)
plt.suptitle("TPM_Clusters")
g.ax_heatmap.set_xticklabels(['FL_Rep1', 'FL_Rep2', 'HL_Rep1', 'HL_Rep2', 'MB_Rep1', 'MB_Rep2'])
plt.savefig('tpm_clustermap.png')
```

## Extra Credit 2
Load the data into DAVID and output the Functional Annotation Chart.  

## Loading more info to IGV  
After downloading the phyloP track from [UCSC Genome Browser](genome.ucsc.edu) and loading it to IGV, we can answer the following questions.  

**What regions seem to have highest PhyloP scores?**  
Exonic regions seem to correlate with peaks in Phylop scores.  

**Are there any highly conserved regions that are not protein-coding (i.e. in exons)? Hypothesize what those might correspond to.**  
There are some high Phylop scores in intronic regions. I can hypothesize that these are pivotal enhancer or activator regions, without which the gene cannot create enough protein product to complete its function.  

## Zooming in on ZRS  

![region of interest](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/regionOfInterest.PNG)

**Take note of the histone modification and conservation patterns at this locus.**  
FL-H3K4me1  1.2 - 1.5 (on a 0 - 3.5 scale)
FL-H3K2ac 0.4 - 0.6 (on a 0 - 7.5 scale)
HL-H3K4me1 1.3 - 1.6 (on a 0 - 4.4 scale)
HL-H3K2ac 1.5 - 3.5 (on a 0 - 14.5 scale)
ML-H3K4me1 0.3 - 0.5 (on a 0 - 3.5 scale)  
ML-H3K2ac 0.3 - 0.5 (on a 0 - 12.4 scale)
Phylop score consistently high  

The modification and conservation patterns were consistently high or low, depending on the tissue and the type of histone modification.  
**Is it well conserved across species?**  
Yes, the Phylop score is high.  

**Based on the histone modifications, for which tissues does this look like a putative enhancer region?**  

This could be a putative enhancer for the midbrain because the methylation levels are low for only this tissue, and methylation is involved in gene repression. Because we see few histone modifications for the midbrain tissue, that area of the genome is more accessible for further processing.  

## Multiple sequence alignment of enhancer sequences  
```
mafft /home/linux/ieng6/cs185s/public/week4/zrs_sequences.fa > mafft_output.fa
```

**What do the "-" characters mean in the MSA output?**  
The '-' characters represent gaps in the alignment, i.e. where we use a character from the reference and not from the individual sequence.  

```
mview -html full -coloring any mafft_output.fa > mafft_output.html
```

**Do you notice any regions that are conserved in all species except snakes? Take a screenshot of those regions.**  
Yes, I notice a couple regions that are conserved for all species except snakes (the last 4 rows in my screenshots).  

![1](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/zrs_1.PNG)  
![2](https://github.com/cse185-sp18/cse185-week4-aarthivenkat/blob/master/labreport/zrs_2.PNG)  

**Based on what we've learned about enhancer regions, hypothesize why this deletion might lead to a loss of legs in snakes.**   
Perhaps this enhancer is a pivotal region for limb development, such that without this enhancer region, the genes for limb development will not be activated, or will not produce enough protein product to develop limbs.  
