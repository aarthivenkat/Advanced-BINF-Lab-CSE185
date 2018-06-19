# Week 5: What color should my eyes be?  (part 1)
Skills covered: GWAS, plink, plotting

The goal of this week's lab is to determine the genetic basis of eye color.

Today, you'll perform a genome-wide association study (GWAS) for samples with blue vs. brown eyes. Thursday, you'll implement a program to predict eye color from a handful of SNPs. Note: if you have done 23andMe, you can even do Part II on your own genome to predict what color your own eyes should be, and see how close that prediction is to reality! 

## 0. Introduction to plink

For our GWAS, we've gone out and collected DNA samples and recorded eye color for a cohort of samples. We've analyzed our samples using a SNP genotype array and cleaned our dataset to remove duplicated samples and problematic markers. We're ready to conduct our GWAS using the cleaned dataset which can be found in the `public/week5` directory with prefix `lab5_gwas_eyecolor`.

Today we'll primarily be using the [plink](https://www.cog-genomics.org/plink2) package, which is a general toolkit for doing all kinds of operations on genetic datasets. Before we dive into GWAS, let's get familiar with the types of files used by `plink`. Almost every `plink` command takes `--file` as an argument, which gives a prefix to the following files:

* `$PREFIX.ped`: This file contains all the genotype information. There is one row per individual. The columns are described [here](https://www.cog-genomics.org/plink/1.9/formats#ped). There are 2V+6 fields, where V is the number of variants. The first 6 columns contain: sample id, family id, id of father (if known), id of mother (if known), sex code, and phenotype value. Each column after that gives the genotype for each variant in each sample. There are two columns giving the two alleles of each sample. Here, 1=reference allele, 2=alternate allele, and 0=missing genotype (this is a little bit different notation than the examples we had in class).
* `$PREFIX.map`: This file describes the variants whose genotypes are given in the `.ped` file. It has four columns: chromosome, SNP ID, position in centimorgans (or 0 if ignored), and genomic position.
* `$PREFIX.fam`: Same as the first 6 columns of the `.ped` file.

We have one additional file with phenotype information, which could have alternatively been stored in the `.ped` file:

* `$PREFIX.phen`: This file has 3 columns: sample id, family id, and phenotype. Here we have coded "cases" (blue eyes) as 2 and "controls" (brown eyes) as 1.

Note since we are not using any family information, we have just used the sample id as the family id.

`plink` alternatively can use binary files by specifying `--bfile` for commands instead of `--file`. In that case, you would have `.bed`, `.bim`, and `.fam` files. It is usually a good idea to work with compressed data when possible to save space. Here we have used the text files to give us a chance to look at the file format on the command line more easily.

Take a look at our dataset. Record:
* How many samples are in our cohort?
* What is the frequency of brown vs. blue eyes in our samples?
* How many markers are included in our dataset?
Record the commands you used to determine this in your lab notebook. Note the results in the methods section of your lab report.

<blockquote>
**UNIX TIP**: In order to count the number of cases and controls in your dataset, you can use unix command `uniq -c`. When applied to a **sorted list**, this command finds the unique values and number of occurence of each value in that list. In order to generate a list of phenotypes, you can use `cat $file | cut -fN` to isolate the `N`'th column of a file. Use the following template to count the number of cases and controls:
    
    `cat $file | cut -fN | sort | uniq -c`
    
</blockquote>

## 1. Analyzing population structure

It turns out our sample was collected from individuals with different ethnic backgrounds. As we will learn in class on Wednesday, using a heterogenous sample may introduce spurious signals in our GWAS. We'd like to control for population structure in our cohort. For this, we'll use something called Principal Component Analysis (PCA). PCA looks for groups of features (in this case, SNPs) that explain the most variation in our data. You can think of it as clustering our samples based on their ancestry. We'll use the results of this clustering as covariates in our GWAS. If this is confusing, revisit this section after Wednesday's lecture!

We can use `plink` to calculate principle components in our sample:

```
plink \
    --file ${PREFIX} \
    --pca 10 \
    --out ${OUTRPREFIX}
```

This will create an output file `${OUTPREFIX}.eigenvec`, with the value along each principal component for each sample for each of the first 10 principle components. Make a scatter plot of PC1 vs. PC2. You can make the plot using your favorite method (e.g. Python, R, or even Excel are all fine). Do you see any distinct clusters of samples? What do you think those clusters correspond to?

## 2. Performing a basic GWAS

Now, we'll use `plink` to perform our GWAS. The following command performs a basic case-control GWAS using logistic regression:
```
plink \
      --file PATH_TO/lab5_gwas_eyecolor \
      --pheno PATH_TO/lab5_gwas_eyecolor.phen \
      --out $OUTPREFIX \
      --logistic \
      --allow-no-sex
```
Note you'll have to fill in absolute paths.

Look at the [`plink` documentation](https://www.cog-genomics.org/plink/1.9/assoc#linear) to learn how to add covariates to our analysis. (Note: the parameter to include the covariates in the analysis, yet remove covariate-specific lines from the output results report.) Run the GWAS twice: once with no covariates and once controlling for population structure using the PCs generated above. Be sure to use a different `$OUTPREFIX` each time so you don't overwrite the original results.

These commands will create an output file `$OUTPREFIX.assoc.logistic`. See the [plink documentation](https://www.cog-genomics.org/plink/1.9/formats#assoc_linear) for a description of each column. There is one row per variant tested. Note, when run with covariates there will also be an additional line of output for each covariate, described under the "Test" column. Use the script below to pull out only the tests for each genotype for downstream analysis (plus keep the header info):

```
cat $OUTPREFIX.assoc.logistic | awk '($5=="ADD" || $0~/CHR/)' > $OUTPREFIX.assoc.logistic.no_covars
```

In each case, how many variants pass genome-wide significance of p<5*10<sup>-8</sup>? Did you get more or fewer significant variants after controlling for covariates?

## 3. Visualizing GWAS results
Now, we'd like to visualize our results. We will use the [`assocplots`](https://github.com/khramts/assocplots) python package for this.  A script for plotting has been provided in `scripts/gwas_plotter.py`. Figure out how to run this script to generate QQ plots and Manhattan plots for each GWAS you performed (with and without covariates). Hint, it is always a good idea to run a script with no arguments to see if there are any hints on how to use it (e.g. `python ./scripts/gwas_plotter.py`)

Include the figures in your lab report. How did the two results differ? Which GWAS do you think is more reliable and why? **Choose that one for the remainder of analyses.**

*Hint*: Take a look at Box 2 in the required reading assignment for more info on how to interpret QQ plots.

## 4. Analyzing significant hits

You will notice some "high rises" in your Manhattan plot. As we discussed in class, it is unlikely that every variant in that spike is independently associated with eye color. It is far more likely that each spike represents a single causal variant which is correlated with a lot of nearby variants that are actually not contributing at all. We can use `plink` to "clump" variants into likely independent signals:

```
plink \
    --file ${PREFIX} \
    --clump ${OUTPREFIX}.assoc.logistic --clump-field P \
    --clump-p1 0.0001 \
    --clump-p2 0.01 \
    --clump-r2 0.5 \
    --clump-kb 250 \
    --out ${OUTPREFIX}
```
Take a look at plink documentation to learn what the different parameters mean.

This command will output a table `$OUTPREFIX.clumped` with a format described [here](https://www.cog-genomics.org/plink/1.9/formats#clumped). How many signals were identified? Report significant hits (meeting genome-wide significance) in a table in your lab report. Use IGV, the UCSC Genome Browser, or another means to determine whether each variant falls within a gene. If it does, include the name of the gene in the table.

Previously, the following SNPs have been associated with eye color:

| chrom | start | rsid | Minor allele |
|----------|----------|-------|------|
| 15 | 28365618 | rs12913832 | A |
| 15 | 28230318 | rs1800407 | T |
| 14 | 92773663 | rs12896399 | G |
| 5 | 33951693 | rs16891982 | C |
| 11 | 89011046 | rs1393350 | A |
| 6 | 396321 | rs12203592 | T |

Do any of these show up in your results? If not, what p-values did you calculate for each one?
Discuss in your lab report why you might not have been able to identify them. We'll use this set of SNPs on Thursday to predict eye color in samples where the answer is unknown.

**That's it for today. Next time, we'll use GWAS hits to *predict* eye color given an individual's genotypes!**
