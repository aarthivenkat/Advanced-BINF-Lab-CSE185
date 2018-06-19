# Week 5: What color should my eyes be?  (part 2)
Skills covered: Trait prediction, scripting, and plotting

Today, we'll use known SNPs associated with eye color to try and predict the eye color of a set of samples.

We will use a set of SNPs previously published in the [IrisPlex paper](https://www.ncbi.nlm.nih.gov/pubmed/20457092) (see the table below). We'll use a slightly more complicated model than a simple case/control model to account for the fact that some people have neither brown nor blue eyes, which we'll denote as "other" here.

## 5. Eye color prediction model and data
The IrisPlex paper has fit the following two models, which can be thought of as performing two separate case control analyses: one of blue vs. brown, and the second of other vs. brown.

<a href="https://www.codecogs.com/eqnedit.php?latex=\ln(p_{blue}/p_{brown})&space;=&space;\alpha_1&space;&plus;&space;\sum_{k}\beta_{1,&space;k}X_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ln(p_{blue}/p_{brown})&space;=&space;\alpha_1&space;&plus;&space;\sum_{k}\beta_{1,&space;k}X_k" title="\ln(p_{blue}/p_{brown}) = \alpha_1 + \sum_{k}\beta_{1, k}X_k" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\ln(p_{other}/p_{brown})&space;=&space;\alpha_2&space;&plus;&space;\sum_{k}\beta_{2,&space;k}X_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ln(p_{other}/p_{brown})&space;=&space;\alpha_2&space;&plus;&space;\sum_{k}\beta_{2,&space;k}X_k" title="\ln(p_{other}/p_{brown}) = \alpha_2 + \sum_{k}\beta_{2, k}X_k" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha_i" title="\alpha_i" /></a>
is the intercept term for model
<a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i" title="i" /></a>, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\beta_{i,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_{i,k}" title="\beta_{i,k}" /></a>
is the effect size for model 
<a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i" title="i" /></a> 
at the 
<a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>
th SNP, and 
<a href="https://www.codecogs.com/eqnedit.php?latex=X_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_k" title="X_k" /></a> 
is the number of minor alleles at an individual's genotype for the <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>th SNP (0, 1, or 2). 

We can rearrange this model to predict the probability of each eye color for a given individual:

<a href="https://www.codecogs.com/eqnedit.php?latex=p_{blue}&space;=&space;\frac{e^{\alpha_1&space;&plus;&space;\sum_k&space;\beta_{1,k}X_k}}{1&plus;e^{\alpha_1&space;&plus;&space;\sum_k&space;\beta_{1,k}X_k}&plus;e^{\alpha_2&space;&plus;&space;\sum_k&space;\beta_{2,k}X_k}}&space;\\\\&space;p_{other}&space;=&space;\frac{e^{\alpha_2&space;&plus;&space;\sum_k&space;\beta_{2,k}X_k}}{1&plus;e^{\alpha_1&space;&plus;&space;\sum_k&space;\beta_{1,k}X_k}&plus;e^{\alpha_2&space;&plus;&space;\sum_k&space;\beta_{2,k}X_k}}&space;\\\\&space;p_{brown}&space;=&space;1-p_{blue}-p_{other}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_{blue}&space;=&space;\frac{e^{\alpha_1&space;&plus;&space;\sum_k&space;\beta_{1,k}X_k}}{1&plus;e^{\alpha_1&space;&plus;&space;\sum_k&space;\beta_{1,k}X_k}&plus;e^{\alpha_2&space;&plus;&space;\sum_k&space;\beta_{2,k}X_k}}&space;\\\\&space;p_{other}&space;=&space;\frac{e^{\alpha_2&space;&plus;&space;\sum_k&space;\beta_{2,k}X_k}}{1&plus;e^{\alpha_1&space;&plus;&space;\sum_k&space;\beta_{1,k}X_k}&plus;e^{\alpha_2&space;&plus;&space;\sum_k&space;\beta_{2,k}X_k}}&space;\\\\&space;p_{brown}&space;=&space;1-p_{blue}-p_{other}" title="p_{blue} = \frac{e^{\alpha_1 + \sum_k \beta_{1,k}X_k}}{1+e^{\alpha_1 + \sum_k \beta_{1,k}X_k}+e^{\alpha_2 + \sum_k \beta_{2,k}X_k}} \\\\ p_{other} = \frac{e^{\alpha_2 + \sum_k \beta_{2,k}X_k}}{1+e^{\alpha_1 + \sum_k \beta_{1,k}X_k}+e^{\alpha_2 + \sum_k \beta_{2,k}X_k}} \\\\ p_{brown} = 1-p_{blue}-p_{other}" /></a>

Below is a reproduced table of the model parameters from the 6 predictive SNPs, converted to hg19 coordinates:

| chrom | start | rsid | Minor allele |
|----------|----------|-------|------|
| 15 | 28365618 | rs12913832 | A |
| 15 | 28230318 | rs1800407 | T |
| 14 | 92773663 | rs12896399 | G |
| 5 | 33951693 | rs16891982 | C |
| 11 | 89011046 | rs1393350 | A |
| 6 | 396321 | rs12203592 | T |

For example, see the following example calculation:
https://docs.google.com/spreadsheets/d/1aP4OQdNsBj7gN5v_Hb40VFV3zgToDPMvGTFQYD6yHwg/edit?usp=sharing

<blockquote>
 **NOT UNIX TIP!** You can do *a lot* with spreadsheets! Although some bioinformaticians might say this is not cool, I think Excel is a great resource especially for working out how to do certain calculations. Please take a chance to look around at this spreadsheet and see how the different fields were calculated. You can actually do all of the rest of this assignment in Excel with some clever rearranging of the fields in the spreadsheet above!
</blockquote>

For this part of the lab, we'll be working with a smaller set of samples independent from our original GWAS. A VCF file containing these variants for our samples can be found in the `public/week5` directory:

```
lab5_pred_eyecolor.vcf.gz
```

We've seen VCF files in the first couple weeks. But let's take a second to remind ourselves what's going on here and point out a couple of important fields. Recall:

* Lines at the top beginning with `#` give header info. The last line of the header gives headers for the different columns describing each variant (columns 1-9), plus one column for each sample in our dataset (columns 10+).
* Each row is for a single variant. Columns 10+ give the genotype (plus other info) for each sample at each position. 
* Note the `ID` columns, which corresponds to the same [dbSNP](https://www.ncbi.nlm.nih.gov/projects/SNP/) rsids given in the table above. Also note the `REF` and `ALT` columns giving the two different alleles at each locus.

Use UNIX commands we've learned in class to determine (1) the number of variants and (2) the number of samples contained in this VCF file. Record your commands in your lab notebook and report the answers in the methods section of your lab report.

Before we move on to predicting eye color, let's convert this VCF file into a format that will be easier for us to process. Use the following commands to create a file with one row per sample and one column per SNP:

```
# Get a header row listing the samples
bcftools query -l lab5_pred_eyecolor.vcf.gz | datamash transpose | awk '{print "ID\t"$0"\t"}' > lab5_pred_eyecolor.tab
# Extract the ID and sample genotypes for each variant
bcftools query -f "%ID\t[%TGT\t]\n" final/lab5_pred_eyecolor.vcf.gz | sed 's/|//g' >> lab5_pred_eyecolor.tab
# Transpose the file
cat lab5_pred_eyecolor.tab | datamash transpose > lab5_pred_eyecolor_transpose.tab
```

Let's break down these commands. You should also run these one at a time to get a good idea of what each command is doing.

1. The first command gets a header for our file, consisting of "ID" and the name of each sample. There are two really useful tools here:

* `bcftools query` can be used to wrangle a VCF file into pretty much any format you want. Here, we used `-l` to simply list the samples present in the VCF file. Using the `-f` option (below), we can output the fields in any way we want.
* `datamash transpose` will *transpose* a file. That is, it makes the rows into columns and columns into rows. So the input here was a list of one sample per line. Transposing that gives us a single row with one sample per column.

The final `awk` command just appends "ID" as the first column. Note, in awk `$0` refers to the entire line that was input.

2. The second command extracts the genotypes for each sample at each variant.

* `bcftools query` is now being used to output genotype info in a specified format. We used `-f` to tell it how we want things output. `%<fieldname>` refers to a VCF field. So `%ID` means output the `ID` column for each line. Anything in brackets (`[]`) is output *per sample*. So `[%TGT\t]` means to print the genotype of each sample separated by tabs. Read more about `bcftools query` [here](https://samtools.github.io/bcftools/bcftools-man.html#query).
* `sed`, like `awk`, is very useful for manipulating text on the command line. Here, we used `sed` to "find and replace" using the syntax `sed 's/<find>/<replace>/g'` to remove all the `|` characters.

3. We used `datamash transpose` again to get a file with one row per sample and one column per SNP.

Note, these SNPs are sorted by genomic coordinate, so are not in the same order as the SNPs in the table above or in the example spreadsheet! You might want to use `awk` to rearrange the columns before moving forward.

## 6. Eye color prediction

Predict the eye color of each sample in our dataset using the model given above. We'll let you decide the best way to do this. For example, you could write a python script. Alternatively, it would be pretty easy to modify the example worked out in the spreadsheet above to do all of this in a spreadsheet with one sample per row. Whichever method you choose, be sure to note in your lab notebook what you did and upload any additional scripts (or spreadsheets) to your Github repository.

For each sample, you should report p(brown), p(blue), and p(other), in addition to a prediction based on what eye color has highest probability for that sample. Include a table of results as a *supplementary file* that you reference from your lab report (do not include the table directly in the report!).

## 7. Comparison of predictions by population

It turns out our sample has a mixture of individuals from two population groups. One (denoted CEU) is from Northern Europe, whereas the second (denoted TSI) is from Southern Europe.

Population labels for each sample are given in `public/week5/ceu_tsi_population_labels.txt`.

Calculate the mean probability of blue, brown, or other colored eyes for each population, CEU and TSI. Which group is more likely to have blue eyes? Does this match with what is known about eye color frequencies in those populations? (Note, not all of these samples has genotyopes, so don't worry if some data is missing).

## 8. Variant interpretation

Now, we'd like to understand a little more about what these variants are doing. Look up these six SNPs using either IGV or the UCSC genome browser (be sure to look at human reference genome build hg19). Where do these SNPs fall? Are they in protein coding regions? Do you have any hypotheses about how these regions might affect eye color? You may also try loading additional tracks to these browsers, such as histone modifications, to see if these SNPs overlap predicted regulatory regions. Some of these can be directly accessed on IGV by going to "File->Load from Encdoe". Keep in mind, that often the SNPs identified by GWAS are not themselves causal variants but might lie nearby truly causal mutations. Discuss your hypotheses in the lab report.

## 9. **(Optional, but fun!) Extra credit (1 pt) **

For extra credit, we can use what we learned about eye color to predict what eye color a child will have based on following parent genotypes:

| chrom | start | rsid | Minor allele | Mother Genotype | Father Genotype |
|----------|----------|-------|------|--------|--------|
| 15 | 28365618 | rs12913832 | A | GG | AG |
| 15 | 28230318 | rs1800407 | T | CC | CC |
| 14 | 92773663 | rs12896399 | G | GG | GT |
| 5 | 33951693 | rs16891982 | C | GG | GG |
| 11 | 89011046 | rs1393350 | A |GA | GA |
| 6 | 396321 | rs12203592 | T | CC | CC |

Based on the fact that each child will inherit one allele from mom and one from dad, figure out all possible sets of SNP genotypes for a child. What is the probability the child will have blue eyes?

## 10. For your lab report
Similar to last week, the lab report document contains specific prompts that should be answered.

**Congrats, you made it to the end of lab 5! This will be our last full lab report before the final project. Next week's lab will just require answering a series of questions to give you a break before we head into final project mode.**
