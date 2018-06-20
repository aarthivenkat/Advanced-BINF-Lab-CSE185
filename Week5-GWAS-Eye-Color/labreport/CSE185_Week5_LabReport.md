# Eye color SNP associations and probalistic modeling through GWAS and literature  
##### Aarthi Venkat | 05 May 2018  

## Abstract

We performed a GWAS for blue versus brown eyes, performing Principal Component Analysis, visualizing Manhattan/QQ Plots, and controlling for population structure. We could only implicate one of six known SNPs associated with eye color, so we turned to known SNPs for downstream analysis. Using a published model, we predicted eye color for 2504 samples, comparing probabilities between two populations and determining that Northern Europeans are more likely to have blue eyes. Finally, we interpreted each variant through studying gene function, SNP position, and histone modifications. This paper is critical in finding phenotype predictors through population genetics and predictive modeling.  

## Introduction

According to [McCarthy et al](https://drive.google.com/file/d/1_PMXVFLCkT0680A4hZ7nXbZkr8QId61x/view)<sup>[1]</sup>, through the advent of large-scale, high-density GWA studies, there has been a rapid increase in the number of SNPs implicated in predisposition for disease or phenotype. Thus, this week we attempt to leverage advances in GWAS, as well as population genomics tool [PLINK](http://zzz.bwh.harvard.edu/plink/)<sup>[2]</sup>, to implicate SNPs associated with eye color - particularly blue versus brown eyes. However, we were only able to associate two SNP clumps, only one of which included on the list of six known SNPs associated with eye color. Thus, we work with the six known SNPs for the remainder of the analysis, predicting probability of eye color through an existing model, and comparing eye colors between Northern and Sourthern Europeans. Finally, we trace these SNPs through a pathway from variant to melanin production to iris color differentiation, an analysis we achieve through case-by-case study of each variant position, gene function, and histone modification at the SNP.  

Ultimately, the question of determining SNPs associated with blue versus brown eye color, as well as predicting eye color and comparing populations, is one situated in a larger picture of finding significant predictors of human phenotypes. By working with population genomic data and utilizing predictive models, we can implicate SNPs as a predictive measure, whether it be for population analysis as in this case, or as a diagnostic or pre-emptive measure, as in the case for breast or prostrate cancer, both of which linked to variants.  

## Methods

### Dataset description

We could analyze the GWAS cohort prior to utilizing Plink to being association studies. There were 261 samples collected, of which 104 had brown eyes and 157 had blue. 911774 markers were genotyped. We know this information through basic UNIX commands. Later, we used a separate cohort for eye color prediction, as our GWAS dataset had weak power and was not able to capture many common SNPs previously associated with eye color. In this dataset, we looked at 2504 samples and 6 SNPs. We used a [custom awk script](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/scripts/unix_custom_scripts) (line 1 in linked file) to identify the number of samples.  

### GWAS

To analyze the GWAS cohort, we turn to [**PLINK v1.90b5.4**](http://zzz.bwh.harvard.edu/plink/)<sup>[2]</sup>. First, we control for population structure by performing Principal Component Analysis with 10 PCs to look for groups of features that could explain variation in the data (in this case, ancestry). From the eigenvec output file, I made a scatter plot of PCA2 vs. PCA1 using [a custom python script](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/scripts/pca1_vs_pca2.py). Then, we perform a logistic regression given a case/control phenotype (blue eyes vs brown eyes), once including covariates and once without including covariates. From the `assoc.logistic` files, we can derive Manhattan plots and QQ Plots using the given [python script](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/scripts/gwas_plotter.py). We interpret from the QQ Plot that the "with covariate" analysis is better due to fewer confounding factors creating falsely significant results, so we move forward using these results. From the Manhattan plot, it is clear that each spike comprises of many smaller, linked signals, so we need to "clump" the variants into independent signals using PLINK. Using the clump field parameter P as well as the following parameters below from the PLINK documentation, I was able to produce a "clumped" file with two significant signals.  

     --clump-p1 0.0001            Significance threshold for index SNPs  
     --clump-p2 0.01              Secondary significance threshold for clumped SNPs  
     --clump-r2 0.50              LD threshold for clumping  
     --clump-kb 250               Physical distance threshold for clumping  

### Eye color prediction  

In order to create a file with one row per sample and one column per SNP, I ran [two custom commands](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/scripts/unix_custom_scripts) (lines 2 and 3) and transposed the file. Then, to get the order of the SNPs in the output to match the given excel spreadsheet, I ran another [awk command](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/scripts/unix_custom_scripts) (line 4). Extending the analysis from the [spreadsheet](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/EyeColorProbabilities.xlsx) (Sheets 1 and 2) the TAs linked, I was able to calculate the probability of each eye color for each sample. Using Excel commands `VLOOKUP` and `AVERAGE`, I could calculate the mean probabilities of blue, brown and other colored eyes for each population in [Sheet 3](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/EyeColorProbabilities.xlsx). To interpret these variants and how they play a role in predicting eye color, I turned to the [**UCSC Genome Browser (build hg19)**](http://genome.ucsc.edu/cgi-bin/hgGateway)<sup>[3]</sup> and overlayed the data with Histone modifications from ChIP-Seq from ENCODE/Broad Institute. Further, for each SNP in or near a gene, I read about the gene's function, particularly with respect to melanin production, using the [**NIH Genetics Home Reference**](https://ghr.nlm.nih.gov/gene)<sup>[4]</sup>.  

## Results  

### GWAS  

After performing Principal Component Analysis, I plotted PCA2 versus PCA1 on a scatter plot in Figure 1. I see three distinct clusters of samples. Those clusters may correspond to features in the underlying dataset that explain variation that may be confounding factors in association tests. Later we are told that ancestry is a differentiating factor, so this is likely a reason behind the clustering.  

#### Figure 1. PCA2 versus PCA1 Scatter Plot  
![pca](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/pca1_vs_pca2.png).

Without controlling for population structure, 24 variants passed genome-wide significance in GWAS. With controlling, 4 variants passed. After clumping, 37 independent signals remained, 2 of which were significant at threshold 5e-8.

Figures 2 and 3 below visualize the QQ plots and Manhattan plots for the GWAS without covariates and with, respectively.  

#### Figure 2. Plots Without Covariates  
##### QQ Plot  
![qq](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/lab5_nocov_qq.png)  
##### Manhattan Plot    
![m](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/lab5_nocov_manhattan.png)  

#### Figure 3. Plots With Covariates   
##### QQ Plot  
![qq](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/lab5_cov_qq.png)  
##### Manhattan Plot    
![m](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/lab5_cov_manhattan.png)  

The QQ Plot differences indicate that the assocation without covariates shows an inflation of assocations. Because the plot is always above the diagonal, we can infer there are confounding factors that are influencing our study of associations. Thus, the GWAS that is more reliable is the one with covariates included in the analysis: the `lab5_cov` files. I used these for the remainder of my analysis.  

After clumping, 2 signals, noted in Table 4, had genome-wide significance.  

#### Table 4. Clumped Significant Signals  
CHR|F|SNP|GENE|BP|P|SP2  
---|---|---|---|---|---|:---  
15|1|rs1129038|HERC2|28356859|2.75e-13|rs12913832(1)  
15|1|rs1470608|OCA2|28288121|3.56e-08|rs749846(1),rs3794604(1),rs4778232(1),rs1448485(1),rs16950821(1),rs8024968(1),rs7177686(1),rs6497253(1),rs7170869(1),rs1375164(1),rs6497254(1)  

### Eye color prediction   
The spreadsheet (Sheet 1) of eye color predictions for each sample can be found [here](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/EyeColorProbabilities.xlsx).  

With the nationality of each sample, I could compute the mean probabilities for each eye color for Northern and Southern Europe, results noted in Table 5. From this analysis, it is evident that CEU individuals are more likely to have blue eyes, which matches our knowledge of eye color frequencies based on lecture slides.  


#### Table 5. Mean Eye Color Probabilities for CEU versus TSI  
Population|P(blue)|P(green)|P(brown)  
---|---|---|---  
CEU|0.631|0.130|0.239   
TSI|0.250|0.167|0.583  

For each SNP, I found the gene it is located in or near, and the histone modifications through the UCSC Genome Browser. For example, for rs1393350, I evaluate the SNP using the visualization in Figure 6. I interpret that rs1393350 is in the intronic region of the gene TYR. According to the [NIH Genetics Home Reference](https://ghr.nlm.nih.gov/gene/TYR), this gene is responsible for producing the enzyme tyrosinase, which converts tyrosine to dopaquinone in the first step of melanin production in melanocytes. Therefore, this SNP is likely associated with a mutation in producing tyrosinase, which may in turn change the levels of melanin production which could differentiate iris color. Because there are low levels of histone modifications (as evidenced by Figure 6), it is difficult to implicate this region as an enhancer with certainty, but I can posit that the region is linked to an enhancer or plays a similar role in the levels of melanin production.  

#### Figure 6. UCSC Genome Browser - rs1393350  
![ucsc](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/ucsc.PNG)  

## Extra Credit  

Based on the parental genotypes, we can determine the probabilities of the child's genotype for each SNP, noted in Table 7.

#### Table 7. Probabilities of Child Genotype for each SNP  

chr|start|rsid|child genotype|probability of each genotype  
---|---|---|---|---  
15|28365618|rs12913832|GA,GG|50%,50%  
15|28230318|rs1800407|CC|100%  
14|92773663|rs12896399|GG,GT|50%,50%   
5|33951693|rs16891982|GG|100%    
11|89011046|rs1393350|GA,GG,AA|50%,25%,25%  
6|396321|rs12203592|CC|100%  

Therefore, there are 16 possible sets of SNP genotypes the child may have (recognizing that for rs1393350, we can split GA into GA and AG, each with 25% probability). Each of the 16 sets have equal probability. These sets are elucidated in the spreadsheet linked ([Sheets 1 and 2](https://github.com/aarthivenkat/Advanced-BINF-Lab-CSE185/blob/master/Week5-GWAS-Eye-Color/labreport/EyeColorProbabilities.xlsx)). We can calculate the probability of the 16 options and average the possibilities for each eye color probablity to result in the following:  

**P(blue) = 0.539843793**  
**P(green) = 0.150725997**  
**P(brown) = 0.30943021**  

Thus, to answer the question in the instructions - **What is the probability the child will have blue eyes?** - 54% likelihood. 

## Discussion  

### Controlling for Covariates

I got fewer significant variants after controlling for covariates, yet these are likely better datapoints, as without controlling, I procured many results that were only significant due to the confounding factors. Through visualization of the QQ plots, it was evident that, without controlling, there was an inflation of associations. The plot was always above the diagonal, and because it is unlikely that the entire GWAS was significant, we can infer that there are confounding factors that are influencing our study of associations. Thus, controlling for population structure, i.e. working with the GWAS with covariates included in the analysis, helped mitigate the effects of these confounding factors.  

### Issues with my GWAS in Identifying Known SNPs  

I only identified one of the known SNPs contributing to eye color - rs12913832. The six given SNPs are associated with variants in eye color, but our dataset is small and only had blue and brown eye colors. Perhaps the SNPs given involve eye color differences between a spectrum of colors, or the differences between blue and green, and our dataset could not capture these SNPs. Further, because our QQ plot for the "with covariate" analysis was plotted largely below the linear trend, we know our data has low power and thus is susceptible to miss some significant results.  

### Hypotheses on links of SNPs to eye color  

**rs12913832** is in the intronic (not protein-coding) region of HERC2 and shows elevated levels of H3K36me3 modifications in nearly all cell lines from ENCODE ChIP-Seq. According to the [NIH Genetics Home Reference](https://ghr.nlm.nih.gov/gene/HERC2), this gene's function is protein ubiquitination, which is a process in protein regulation. Further, genetic variations in this gene are associated with skin/hair/eye pigmentation variability. Thus, I can hypothesize that this SNP is either causal or linked with a causal mutation associated with changing protein ubiquitination activity, such that the protein product of OCA2 and other genes HERC2 regulates has differential activity for different SNPs. Because [OCA2 produces a protein responsibile for producing the pigment melanin](https://ghr.nlm.nih.gov/gene/OCA2), perhaps when there is decreased activity of OCA2, there is a direct association with lighter eyes.  

**rs1800407** is in an exonic (protein-coding) region OCA2. Thus, perhaps the SNP in the DNA encodes (or is near a mutation that encodes) for a different amino acid in the peptide sequence, which can change the conformation of the protein producing the pigment for eye color. Thus could result in a predictively different iris color phenotype.  

**rs12896399** is 16,489bp upstream (not protein-coding) of the gene SLC24A4, which plays a role in potassium-dependent calcium transport. This gene is not as strongly implicated with eye color in terms of function, although we could potentially link calcium transport to transport of the proteins for the production of melanin. Another concept to think about is the possibility that the SNP is at a position close to HERC2 and OCA2 in 3D space, which we could confirm using Hi-C chromatin conformation capture data. If this is the case, perhaps the SNP, whether causal or near a causal variant, could change protein levels of the OCA2 gene through 3D enhancer interactions, which would in turn change eye color.  

**rs16891982** is in the exonic (protein-coding) region of the gene SLC45A2, which, according to the [NIH Genetics Home Reference](https://ghr.nlm.nih.gov/gene/SLC45A2), provides instructions for making a protein that is located in melanocytes, which produce pigments such as that for eye color. The exact function of the gene is likely in transport due to being situated in the class of solute carriers like SLC24A4. Thus, perhaps this SNP encodes (or is near a mutation that encodes) for a different amino acid in the transport protein, which can change the conformation of the protein such that it interrupts the melanin pathway and differentiates pigment levels.  

**rs1393350** is in the intronic (not protein-coding) region of the gene TYR with low levels of histone modifications. [This gene is responsible for producing the enzyme tyrosinase](https://ghr.nlm.nih.gov/gene/TYR), which converts tyrosine to dopaquinone in the first step of melanin production in melanocytes. Therefore, this SNP is likely associated with a mutation in producing tyrosinase, which may in turn change the levels of melanin production which could differentiate iris color.  

**rs12203592** is in the intronic (not protein-coding) region of the gene IRF4 with high levels of various methylation histone modifications (e.g. H3K4me3, H3K27me3). Because of these histone modifications, we can infer that this region is an enhancer region for either IRF4 or genes with which IRF4 works in conjunction. The function of IRF4 is [regulation of interferons in response to infection by virus](https://ghr.nlm.nih.gov/gene/IRF4), which is not as tightly linked to melanin as many of the other SNPs. Perhaps this SNP is linked to a gene that is part of the melanin pathway.  

## Citations

1. McCarthy, M. I. et al. Genome-wide association studies for complex traits: consensus, uncertainty and challenges. Nature Rev. Genet. 9, 356–369 (2008). An informative overview of key issues in the field of GWA studies.  
2. Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.  
3. Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, Haussler D. The human genome browser at UCSC. Genome Res. 2002 Jun;12(6):996-1006.  
4. National Library of Medicine (US). Genetics Home Reference [Internet]. Bethesda (MD): The Library; 2013 Sep 16 [cited 2018 May 05]. Available from: https://ghr.nlm.nih.gov/.  

