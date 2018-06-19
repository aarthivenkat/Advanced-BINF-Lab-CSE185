# CSE 185 Lab Notebook - Week 5  
#### Aarthi Venkat  
#### Date: May 1 2018  

## 0. Introduction to Plink  
**How many samples are in our cohort?**  
`wc -l lab5_gwas_eyecolor.ped`  
261  

**What is the frequency of brown vs. blue eyes in our samples?**  
`cat lab5_gwas_eyecolor.phen | cut -f 3 | sort| uniq -c`  
count color (1=brown, 2=blue)  
104 1  
157 2  

**How many markers are included in our dataset?**  
`wc -l lab5_gwas_eyecolor.map`  
911774  

## 1. Analyzing population structure  
In the public week 5 directory,  
`plink --file lab5_gwas_eyecolor --pca 10 --out ../../avenkat/week5/lab5_gwas_eyecolor`  

From output file lab5_gwas_eyecolor.eigenvec, I made a scatter plot of PCA2 vs. PCA1 using [a custom python script](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/scripts/pca1_vs_pca2.py).  
![pca](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/pca1_vs_pca2.png)  

**Do you see any distinct clusters of samples? What do you think those clusters correspond to?**  
I see three distinct clusters of samples. Those clusters may correspond to features in the underlying dataset that explain variation that may be confounding factors in association tests.  

## 2. Performing a basic GWAS  
**Without covariates**  
`plink --file /home/linux/ieng6/cs185s/public/week5/lab5_gwas_eyecolor --pheno /home/linux/ieng6/cs185s/public/week5/lab5_gwas_eyecolor.phen --out /home/linux/ieng6/cs185s/avenkat/week5/lab5_nocov --logistic --allow-no-sex`  

**With covariates**  
`plink --file /home/linux/ieng6/cs185s/public/week5/lab5_gwas_eyecolor --pheno /home/linux/ieng6/cs185s/public/week5/lab5_gwas_eyecolor.phen --out /home/linux/ieng6/cs185s/avenkat/week5/lab5_cov --logistic --allow-no-sex --covar lab5_gwas_eyecolor.eigenvec --hide-covar`  

Because I used the --hide-covar flag, I don't need to run this command.  
`cat $OUTPREFIX.assoc.logistic | awk '($5=="ADD" || $0~/CHR/)' > $OUTPREFIX.assoc.logistic.no_covars`  

**In each case, how many variants pass genome-wide significance of p<5x10-8?**  
**Without covariates**  
`cat lab5_nocov.assoc.logistic | awk '($9<0.00000005)' | wc -l`  
24  
**With covariates**  
`cat lab5_cov.assoc.logistic | awk '($9<0.00000005)' | wc -l`  
4  

**Did you get more or fewer significant variants after controlling for covariates?**
I got fewer significant variants after controlling for covariates, yet these are likely better datapoints as I controlled for false signficance (confounding factors). 

## 3. Visualizing GWAS Results  
**Without covariates**  
`python ./scripts/gwas_plotter.py lab5_nocov.assoc.logistic lab5_nocov`  
##### QQ Plot  
![qq](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/lab5_nocov_qq.png)  
##### Manhattan Plot    
![m](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/lab5_nocov_manhattan.png)  

**With covariates**  
`python ./scripts/gwas_plotter.py lab5_cov.assoc.logistic lab5_cov`  
##### QQ Plot  
![qq](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/lab5_cov_qq.png)  
##### Manhattan Plot    
![m](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/lab5_cov_manhattan.png)  

**How did the two results differ? Which GWAS do you think is more reliable and why?**  
The QQ Plot differences indicate that the assocation without covariates shows an inflation of assocations. Because the plot is always above the diagonal, we can infer there are confounding factors that are influencing our study of associations. Thus, the GWAS that is more reliable is the one with covariates included in the analysis: the `lab5_cov` files. I will be using these for the remainder of my analysis.  

## 4. Analyzing significant hits  
**Clumping significant hits**  
`plink --file lab5_gwas_eyecolor --clump lab5_cov.assoc.logistic --clump-field P --clump-p1 0.0001 --clump-p2 0.01 --clump-r2 0.5 --clump-kb 250 --out lab5_covclump`  
**How many signals were identified?**  
Output:`--clump: 37 clumps formed from 76 top variants.` 37 signals.  

**Report significant hits (meeting genome-wide significance) in a table in your lab report.**  

`cat lab5_covclump.clumped | awk '($5<0.00000005)'`  

CHR|F|SNP|BP|P|SP2
---|---|---|---|---|:---
15|1|rs1129038|28356859|2.75e-13|1|rs12913832(1)  
15|1|rs1470608|28288121|3.56e-08|rs749846(1),rs3794604(1),rs4778232(1),rs1448485(1),rs16950821(1),rs8024968(1),rs7177686(1),rs6497253(1),rs7170869(1),rs1375164(1),rs6497254(1)  

**Determine whether each variant falls within a gene.**  
Using [**dbSNP**](https://www.ncbi.nlm.nih.gov/snp)  

SNP|gene  
---|---  
rs1129038| HERC2  
rs12913832| HERC2  
rs1470608| OCA2  
rs749846| OCA2  
rs3794604| OCA2  
rs4778232| OCA2  
rs1448485| OCA2  
rs16950821| OCA2  
rs8024968| OCA2  
rs7177686| OCA2  
rs6497253| OCA2  
rs7170869| OCA2  
rs1375164| OCA2  
rs6497254| OCA2  

**Do any of these show up in your results? If not, what p-values did you calculate for each one?**  

chr|start|rsid|Minor allele|Show up in results?|P-val  
---|---|---|---|---|---  
15|28365618|rs12913832| A | YES | 2.75e-13  
15|28230318|rs1800407| T | NO | NOT IN ANY RESULTS  
14|92773663|rs12896399| G | NO | NOT IN ANY RESULTS  
5|33951693|rs16891982| C | NO | NOT IN ANY RESULTS  
11|89011046|rs1393350| A | NO | NOT IN ANY RESULTS   
6|396321|rs12203592| T | NO | NOT IN ANY RESULTS  

**Discuss in your lab report why you might not have been able to identify them.**  
The given SNPs are associated with variants in eye color, but our dataset is small and only had blue and brown eeye colors. Perhaps the SNPs given involve eye color differences between green and a different color. Further, because our QQ plot for the "with covariate" analysis was plotted below the linear trend, we know our data has low power and thus is susceptible to miss some significant results.  

## 5. Eye color prediction model and data  

**Use UNIX commands we've learned in class to determine (1) the number of variants and (2) the number of samples contained in this VCF file**  
1. `zcat *vcf.gz | grep -v "^#" | wc -l` **6 variants**  
2. `zcat lab5_pred_eyecolor.vcf.gz | grep -v "^#" | awk '{print $8}' | awk -F '[ ;]' '{print $11}' | grep "^NS"` **2504 samples**  

**Commands to create a file with one row per sample and one column per SNP**  
```
bcftools query -l lab5_pred_eyecolor.vcf.gz | datamash transpose | awk '{print "ID\t"$0"\t"}' > lab5_pred_eyecolor.tab
bcftools query -f "%ID\t[%TGT\t]\n" final/lab5_pred_eyecolor.vcf.gz | sed 's/|//g' >> lab5_pred_eyecolor.tab
cat lab5_pred_eyecolor.tab | datamash transpose > lab5_pred_eyecolor_transpose.tab
```
To move around the columns to match the spreadsheet,  
`cat lab5_pred_eyecolor_transpose.tab | awk '{print $1, $5, $4, $3, $6, $2, $7}' > temp_lab5_pred_eyecolor_transpose.tab`
`mv temp_lab5_pred_eyecolor_transpose.tab lab5_pred_eyecolor_transpose.tab`  

## 6. Eye color prediction  
Extending the analysis from the spreadhsheet the TAs linked, I was able to calculate the probability of each eye color for each sample. The spreadsheet (Sheets 1 and 2) elucidating the model and outputting the probabilities is linked [here](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/EyeColorProbabilities.xlsx).  

## 7. Comparison of predictions by population  
Using Excel commands `VLOOKUP` and `AVERAGE`, I could calculate the mean probabilities. See my work in [Sheet 3 here](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/EyeColorProbabilities.xlsx).  

**Calculate the mean probability of blue, brown, or other colored eyes for each population, CEU and TSI**  

Population|P(blue)|P(green)|P(brown)  
---|---|---|---  
CEU|0.631|0.130|0.239   
TSI|0.250|0.167|0.583  

**Which group is more likely to have blue eyes? Does this match with what is known about eye color frequencies in those populations?**  CEU is more likely to have blue eyes, which mateches our knowledge of eye color frequencies based on lecture slides (as well as general information about Northern versus Southern Europeans).  

## 8. Variant interpretation  

**Where do these SNPs fall? Are they in protein coding regions? Do you have any hypotheses about how these regions might affect eye color?**  

Using the [UCSC Genome Browser (build hg19)](http://genome.ucsc.edu/cgi-bin/hgGateway) and overlaying the data with Histone modifications from ChIP-Seq from ENCODE/Broad Institute,  

**rs12913832** is in the intronic (not-protein coding) region of HERC2 and shows elevated levels of H3K36me3 modifications in nearly all cell lines from ENCODE ChIP-Seq. According to the [NIH Genetics Home Reference](https://ghr.nlm.nih.gov/gene/HERC2), this gene's function is protein ubiquitination, which is a process in protein regulation. Further, genetic variations in this gene are associated with skin/hair/eye pigmentation variability. Thus, I can hypothesize that this SNP is either causal or linked with a causal mutation associated with changing protein ubiquitination activity, such that the protein product of OCA2 and other genes HERC2 regulates has different activity. Because [OCA2 produces a protein responsibile for producing the pigment melanin](https://ghr.nlm.nih.gov/gene/OCA2), perhaps when there is decreased activity of OCA2, there is a direct association with blue eyes.  

**rs1800407** is in a exonic (protein-coding) region OCA2. Thus, perhaps the SNP in the DNA encodes (or is near a mutation that encodes) for a different amino acid in the peptide sequence, which can change the conformation of the protein producing the pigment for eye color.  

**rs12896399** is 16,489bp upstream of the gene SLC24A4, which plays a role in potassium-dependent calcium transport. This gene is on a different chromosome than HERC2 and OCA2, and is not as strongly implicated with eye color in terms of function. However, perhaps the SNP is at a position close to HERC2 and OCA2 in 3D space, which we could confirm using Hi-C chromatin conformation capture data. If this is the case, perhaps the SNP, whether causal or near a causal variant, could change protein conformation of the OCA2 gene through 3D interactions, which would in term change eye color. I did not see any clear indication of histone modifications from the ChIP-Seq tracks.  

**rs16891982** is in the exonic region of the gene SLC45A2, which, according to the [NIH Genetics Home Reference](https://ghr.nlm.nih.gov/gene/SLC45A2), provides instructions for making a protein that is located in melanocytes, which produce pigments such as that for eye color. The exact function of the gene is likely in transport due to being situated in the class of solute carriers like SLC24A4. Thus, perhaps this SNP encodes (or is near a mutation that encodes) for a different amino acid in the transport protein, which can change the conformation of the protein such that it doesn't complete its function in the melanin pathway to the same extent, which would result in less pigmentation in the iris (e.g. blue eyes).  

**rs1393350** is in the intronic region of the gene TYR with low levels of histone modifications. [This gene is responsible for producing the enzyme tyrosinase](https://ghr.nlm.nih.gov/gene/TYR), which converts tyrosine to dopaquinone in the first step of melanin production in melanocytes. Therefore, this SNP is likely associated with a mutation in producing tyrosinase, which may in turn change the levels of melanin production which could differentiate iris color.  

**rs12203592** is in the intronic region of the gene IRF4 with high levels of various methylation histone modifications (e.g. H3K4me3, H3K27me3). Because of these histone modifications, we can infer that this region is an enhancer region for either IRF4 or genes with which IRF4 works in conjunction. The function of IRF4 is [regulation of interferons in response to infection by virus](https://ghr.nlm.nih.gov/gene/IRF4), which is not as tightly linked to melanin as many of the other SNPs. Perhaps this SNP is linked to a gene that is part of the melanin pathway.  

## 9. Extra Credit  

chr|start|rsid|child genotype|probability of each genotype  
---|---|---|---|---  
15|28365618|rs12913832|GA,GG|50%,50%  
15|28230318|rs1800407|CC|100%  
14|92773663|rs12896399|GG,GT|50%,50%   
5|33951693|rs16891982|GG|100%    
11|89011046|rs1393350|GA,GG,AA|50%,25%,25%  
6|396321|rs12203592|CC|100%  

Therefore, there are 16 possible sets of SNP genotypes the child may have (recognizing that for rs1393350, we can split GA into GA and AG, each with 25% probability). Each of the 16 sets have equal probability. These sets are elucidated in the spreadsheet linked ([Sheets 1 and 2](https://github.com/cse185-sp18/cse185-week5-aarthivenkat/blob/master/labreport/EyeColorProbabilities.xlsx)). We can calculate the probability of the 16 options and average the possibilities for each eye color probablity to result in the following:  

**P(blue) = 0.539843793**  
**P(green) = 0.150725997**  
**P(brown) = 0.30943021**  

Thus, to answer the question in the instructions - **What is the probability the child will have blue eyes?** - 54% likelihood.  
