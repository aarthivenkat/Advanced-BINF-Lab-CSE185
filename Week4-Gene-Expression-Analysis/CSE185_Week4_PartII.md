# Week 4: Why don't snakes have legs? (Part 2)
Skills covered: Data visualization and genome browsers, accessing public databases, epigenomics, multiple sequence alignment

For the first half of the lab today, we'll continue where we left off analyzing expression from different tissues of a developing mouse to identify differentially expressed genes.

We'll then zoom in on one differentially expressed gene, "Sonic Hedgehog" [Shh](https://en.wikipedia.org/wiki/Sonic_hedgehog). This gene is famous for its role in early development. Depending on its concentration in different parts of a developing embryo, it can determine whether cells turn into different structures such as digits, limbs, or brain tissue.

Shh expression is controlled by a neaby enhancer known as the ZRS (Zone of polarizing activity regulatory sequence), which is specific to limb tissues. Intriguingly, mutations in this enhancer in humans have been shown to cause limb malformations. In the second part of the lab today, we'll take a closer look at the ZRS region in mouse, then compare the sequence at this region to other species with and without limbs.

## 6. Differential expression analysis

First we'll use [sleuth](https://pachterlab.github.io/sleuth) to identify differentially expressed genes. `sleuth` is designed to work directly with the output of `kallisto`, which we ran on Tuesday.

We'll need to use R for this. To open the R environment, type:

```
R
```

The following code will run sleuth. Since we don't require previous experience with `R` for the course, we have gone ahead and written these commands for you. However you are encouraged to take a look at the [Sleuth tutorial](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html) to learn more about each of these commands.
```R
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
write.table(sleuth_significant, "~/week4/sleuth_results.tab", sep="\t", quote=FALSE)
```

This will output significant hits to `sleuth_results.tab`. How many significant transcripts are there? Include the results in your lab report.

Take a look at the first couple examples. You'll notice the transcript ID is a big confusing number, e.g. "ENSMUST00000061745.4". To see what gene name that corresponds to, you can go to http://uswest.ensembl.org/Mus_musculus/Info/Index and use the search box in the upper right. For several top hits, find the gene name, navigate to that gene in IGV, and take screenshots to include in your lab report. What are the gene names for the top 10 genes? Be sure to include the gene Shh and the surrounding region in your examples. It should be close to the top of your list. 

**Optional extra credit (1 point)**
Visualize expression of differntially expressed genes as a heatmap (e.g. with transcripts as rows and samples as columns). Cluster the rows by row and column. Do replicates cluster together? Are there clear clusters of up vs. down regulated genes in each tissue?

**Optional extra credit (1 point)**
Perform gene ontology (GO) analysis or gene-set enrichment analysis (GSEA) on the top set of genes. There are many online tools for doing this (e.g. [DAVID](https://david.ncifcrf.gov/), [Panther](http://pantherdb.org/)). What types of biological processes are enriched in differentially expressed genes? Inlude a table in your results section.

## 7. Loading more info to IGV
Launch IGV and load the session you started last Tuesday. You should have already 6 tracks: 2 for each RNA-seq replicate of HL, FL, and MB. Additionally the default Refseq genes should be present at the bottom. We will be focusing today on the Shh region. In the search bar at the top, navigate to region chr5:28,278,817-29,447,265.

We'd like to identify potential regulatory regions for the differentially expressed gene *Shh*. Our labmates have generated some ChIP-sequencing data that will be useful for identifying putative enhancers in this region. We'll start by adding the ChIP-seq datasets and some additional tracks to IGV to help us interpret some features of this region.

In the `public/week4` directory, you'll find (bedGraph)[https://genome.ucsc.edu/FAQ/FAQformat.html#format1.8] files (`*.bedGraph`) for ChIP-sequencing experiments for H3K27ac and H3K4me1 (histone modifications found near enhancer regions). Load the files to IGV. Note these have been restricted to have data in our region of interest, so you if you scroll outside this region you won't see any data.

First take a look at where these marks (H3K27ac and H3K4me1) are falling. Are they near gene regions? Beginning or ends of genes? Other places? Discuss this in the results section of your lab report.

Before moving on, let's add one additional track to IGV about sequence conservation. For this, we'll load the PhyloP track, which gives a per-base pair measure of sequence conservation across species (see Wednesday slides). You can get the PhyloP track from the genome.ucsc.edu (another genome browser! which you can load tracks to similar to IGV.). As part of the genome browser, UCSC hosts many different "tracks" of information for many different genome builds. If you are ever looking for information on any sort of genomic annotation (e.g. gene annotations, conservation, genetic variation, and more), UCSC is a good place to start.

Use the Table Browser (`Tools > Table Browser`) to download the PhyloP track for our region of interest. From the home page, select "Tools->Table Browser". use the dropdown boxes to select the mouse mm10 genome build. Additionally choose the following options:
* Group: Comparative Genomics
* Track: Conservation
* Table: 60 Vert. Cons (phyloP60wayAll)
* Region: position chr5:23414443-36455956
* Filter: edit and select "Limit data output to" 10 million lines otherwise the output will be truncated
* Output format: custom track
* Output file: put a reasonable filename (e.g. `phyloP60wayAll_mm10_Shh_region.wig`)
After you click "get output" make sure "DATA VALUE" is selected on the next page, then click "Get custom track in file" which will download the file to your computer's default download location.

It's worth taking a second to go back to the table browser and see what kinds of info can be downloaded from here. It is a really flexible tool for a huge number of published genome-wide data tracks!

Now go back to IGV and load the PhyloP file. What regions seem to have highest PhyloP scores? Are there any highly conserved regions that are not protein-coding (i.e. in exons)? Hypothesize what those might correspond to. Include a brief description of what you observe in your lab report.

## 8. Zooming in on ZRS

Use IGV to zoom in on region chr5:29,314,718-29,315,770. This region corresponds to the ZRS (Zone of polarizing activity regulatory sequence, also called MFCS1) which is one of the most deeply studied mammalian enhancer sequences known to regulate the *Shh* gene. Take note of the histone modification and conservation patterns at this locus. Is it well conserved across species? Based on the histone modifications, for which tissues does this look like a putative enhancer region? Take a screen shot of this region and include it as well as a description in your lab report.

## 9. Multiple sequence alignment of enhancer sequences

After discussing with your labmates and reading about the ZRS enhancer, you're convinced this is likely an important region involved in limb development. You also recall that the PhyloP scores showed this region is highly conserved. To study how this region differs across species with and without limbs, you look for regions similar to the mouse ZRS in other organisms. Remarkably, you are able to identify similar regions in human, mouse, cow, dolphin (with limbs) and python, rattlesnake, cobra, and boa (snakes without limbs). These are all collected in the file `zrs_sequences.fa` in the `public/week4` folder.

`cat` this fasta file to look at sequences. Each fasta header line gives more info on which genome build and coordinates each region was taken from.

We'd first like to perform multiple sequence alignment (MSA) to see exactly how similar these regions are across species. For that, we'll use the `mafft` tool. Type `mafft --help` at the command line to see usage information for this tool. It should take a fasta file of sequences as input and produce a similar fasta file with the MSA as output. Examine the output file. What do the "-" characters mean?

Now, we will use the Mview tool to visualize the alignment. You can either do this using the web version at https://www.ebi.ac.uk/Tools/msa/mview/ or using the command line `mview` tool. Type `mview --help` to see full usage details. The following command:
```
mview -html full my_msa.fa > my_msa.html
```
will produce a (not very colorful) html file to visualize the MSA. Play around with the options to make the plot nicer and more colorful. Visualize the resulting html file in your web browser.

Do you notice any regions that are conserved in all species except snakes? Take a screenshot of those regions.

You should be able to find at least one region that is deleted in all snakes but conserved across all other species. Based on what we've learned about enhancer regions, hypothesize why this deletion might lead to a loss of legs in snakes. Discuss your hypothesis in your lab report.


## 10. For your lab report

For this week's lab report, there are specific prompts and instructions included in the template document in the `labreports` folder. Each section lists how many points it will be worth, so be sure to complete all the items listed there.
