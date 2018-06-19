# Week 4: Why don't snakes have legs? (Part I)
Skills covered: Genome browsers and data visualization, RNA-seq analysis

You are interested in determining which genes and regulatory regions are most important for controlling the development of limbs (legs and arms), and whether those regions of the genome might explain why some animals are missing limbs (e.g., most snakes). You decide to first tackle this question by analyzing tissues from mouse early development. 

Today, you will use RNA-sequencing data to quantify expression levels and learn how to visualize this data using a genome browser.
On Thursday, you will perform differential expression analysis and use additional epigenomics datasets to further analyze a limb-specific enhancer sequence and trace it's evolution across distantly related species.

You reason that by comparing expression in limb tissues vs. non-limb you can determine genes specific to pathways controlling limb development. You go ahead and collect samples for three tissues from a developing mouse: hind limb (HL), fore limb (FL), and mid-brain (MB). For each tissue, you collect two replicates so you can perform more robust differential expression analysis.

## 1. Data inspection and quality control

First, remember to remove week1-week3 data if you have completed your lab reports for those weeks. 

As in previous weeks, start by making a clone of your repository on `ieng6`:

```
git clone https://github.com/cse185-sp18/cse185-week4-<username>.git week4
```

Data for this week can be found in the `public/week4` directory. Among other things, you can find fastq files for our RNA-sequencing experiments there:

* Forelimb: `FL_Rep1_chr5_*.fq.gz`, `FL_Rep2_chr5_*.fq.gz` 
* Hindlimb: `HL_Rep1_chr5_*.fq.gz`, `HL_Rep2_chr5_*.fq.gz` 
* Midbrain: `MB_Rep1_chr5_*.fq.gz`, `MB_Rep2_chr5_*.fq.gz` 

(Note: these reads have been downsampled from the original experiment to only contain chromosome 5)

First, take a look at the fastq files. **Do not unzip them!** See the UNIX tip below for how to deal with `.gz` files. What read length and how many reads were used for each experiment? Record the answer in your lab notebook and describe in the methods section of your report.

<blockquote>
**UNIX TIP**: Using a compression method like `gzip` or `bgzip` can save tons of space when dealing with huge files. Gzipped files aren't directrly human-readable. However, you can use the `zcat` command to write the contents of the file to standard output. For instance, to see the head of a gzipped file, you can do `zcat file.gz | head`. You can similarly pipe the output of `zcat` to other commands like `wc`.
</blockquote>

Run `fastqc` on the rep1 and rep2 Forelimb, Hindlimb, Midbrain fastq files. You do not need to include the figures in your lab report, although you should keep track of the output `html` files. Comment on anything flagged as problematic in the methods section of your report. See if you can find an explanation in the fastqc help online about whether the flags you see are specific to RNA-seq datasets.

<blockquote>
**UNIX TIP**: If you feel comfortable to do so, then you can create a shell script (.sh) with a variable parameter as the fastq file name and run the shell script in 'screen' or 'nohup' to run the fastqc processing in parallel for all sets of files. See below for more tips on using `screen`. If you do write any of your own scripts for analyses, be sure to include those in the `scripts/` directory of the Github repository for this week and refer to those in your lab notebook.
</blockquote>

## 2. RNA-seq sequence alignment

Since we'd like to visualize where our RNA-seq reads are falling in the genome, we'll need to align the reads to a reference transcriptome. Depending on the alignment method used, this step can either consume tons of memory (e.g. STAR) or take a long time to run (e.g. TopHat). Luckily, you are working with a bioinformatician in the lab who has already gone ahead and aligned the reads to the latest mouse reference for you.

The BAM files containing alignments can be found in the `week4` directory (`*_chr5.bam`).

Unforunately, your labmate forgot to write down some details about how they performed the alignment. Particularly, we'd like to know:

* What program and command did he or she use to align the reads?
* What reference genome were the sequences aligned to?

You wonder if the BAM header files might have some more information about how those files were generated. Type
```
samtools view
```
to learn how to output the header of a BAM file. Take a look at the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) to see a description of different standard header tags (top of page 4). Note the "@PG" tag gives info about the program that was used, and many tools (including the one used here) will use that tag to document the exact command used to generate the BAM file.

Take a look at reads, for instance by doing `samtools view FL_Rep1_chr5.bam`. If you scroll down, you'll notice the CIGAR scores have some extra characters in them we haven't seen before (See week 1 slides for a refresher on CIGAR scores). In the past, we have seen "M" for match, "I" for insertion, and "D" for deletion. Now we see many reads have "N" in the CIGAR scores (e.g. read ID SRR3950230.31710737). In the SAM specification, go to page 6 to read more about CIGAR scores and find out what "N" represents. What biological feature do you think the "N"'s stand for?

<blockquote>
**UNIX TIP**: `less` is really helpful for looking at and scrolling through files. A helpful way to visualize a sam file is to run `samtools view file.bam | less -S`. The `-S` parameter tells the terminal not to wrap lines, and instead allows you to scroll through long lines horizontally. This makes files with long lines much more readable. Another trick: once you're looking at a file using `less`, you can use `ctrl-v` to scroll down more quickly than using the down button.
</blockquote>

## 3. Quantifying gene expressiom

We'll first want to use the RNA-sequencing data to quantify expression of each gene. For this, we'll use a tool called `kallisto` (https://pachterlab.github.io/kallisto/). This is an extremely fast method to quantify transcript abundance. It's main speedup over competing methods is to avoid the alignment step altogether and use a much simpler kmer counting approach based on our old friend from last week: the De Brujn graph! Kallisto will take in our fastq files and output estimated "transcripts per million reads" (TPM) values for each transcript.

To help you get started running `kallisto`, you should find a bash script `run_kallisto.sh` in the Github repository for this week under the `scripts` folder. First take a look at this script to familiarize yourself with what it's doing.

The top of the script sets variables for things each `kallisto` run will use, inclusing the gene annotation (GTF) file and kallisto index (`mm10_kallisto`). The `for` loop loops through each replicate of each experiment and does a separate `kallisto` job.

Type `kallisto quant` to see a description of each option and the syntax for running. Edit the command to use 100 bootstrap samples. Describe any default parameters used in the methods section of your lab report.

You can run the script by:
```
~/week4/scripts/run_kallisto.sh 
```

This may take a while to run (~20 minutes). While you are waiting, move on to part 3 to visualize the RNA-seq data, which can be done independently of the `kallisto` run.

<blockquote>
**UNIX TIP**: To run your bash scripts as executables, you need to change the access premissions of the file to allow execution. To do so, use the following command:
  
  `chmod +x sample.sh `

You can now run your script as an executable: `./script.sh`.
Alternatively, you can run a non-executable bash script with `bash script.sh` command.
</blockquote>

<blockquote>
**UNIX TIP**: It would be efficient to run this script in the 'background' so that you can continue working in your terminal prompt. This can be done with multiple methods such as 'screen' or 'nohup'. When using screen you essentially open a new terminal screen in your current terminal window. To do this type `screen -S kallisto` and you will be in a new screen window. (You can type `pwd` to see where you are.) Now run the script `~/week4/scripts/run_kallisto.sh` and you should see the script running. You can now return to your main terminal window while the script runs in the "background screen", by pressing the keys "control" + "A" + "D" together. Type `screen -ls` to see the screens you have open. In 10 minutes, return to the kallisto screen to check on the progress of your script by typing `screen -r kallisto`.
</blockquote>

## 4. Visualizing data using a genome-browser
Now we'd like to visualize these alignments to help us see which genes might be differentially expressed between our samples. We'll do this statistically in section 6.

For genomic DNA sequences, we previously used `samtools tview` to visualize alignments. This is great if we are looking at genetic variation in one sample, but is less helpful for visualizing *multiple samples* and *read abundances*. Today, we'll introduce a **genome browser** called the [Integrative Genomics Viewer](https://igv.org/), or IGV. On Thursday we'll also use some features of a different genome browser run by UCSC.

To install IGV, go to http://software.broadinstitute.org/software/igv/download and click "Launch with 750MB". Follow the dialog prompts, which should result in an "IGV" icon being created on your Desktop. Double click the IGV icon to launch the program.

After you launch IGV, you'll need to tell it which reference genome to use. In the top left, choose the genome-build you determined was used above. If you're not sure, ask your TA or friend before moving on since nothing will make sense in IGV if you're not using the right genome build!

Take a moment to orient yourself with IGV. It is basically like a Google Maps for genomes! The top gives the names of each chromosome. The bottom track, labeled "Refseq genes" gives the names and coordinates for all annotated genes. Let's choose one to look at. Type "Nanog" in the search box at the top. This will zoom the view in on this gene. Notice how in the gene you can see the exon and intron structure. The little arrows in the introns point to the right, which means this gene is on the forward strand of the reference. Take a look at another gene (e.g. Sox5) to see a gene on the reverse strand. Drag your mouse over different coordinate windows to zoom in further until you can see actual DNA sequence at the bottom. 

Now we'd like to load our sequence alignments. While IGV can directly visualize BAM files, we'll instead look at "counts" files (".tdf" format) which are much smaller. These give read counts per position (i.e. coverage) which can give us an idea of the abundance of each gene. `igvtools` was already used to create `.tdf` files for you (from the entire genome, not just chromosome 5) in the `week4` directory. From IGV, select "File->Load from File" and select the 6 tdf files from the `public/week4` directory to upload. If you are using your own computer, you'll first need to `scp` the files to your local computer before loading them with IGV.

Navigate to a gene. A good one is "chr3:29,939,546-30,023,181" (the gene Mecom). Note that the RNA-seq tracks have very "spiky" coverage. Some regions have tons of reads and others are flat. Note how that compares to the structure of the gene annotated on the bottom. As expected, the "spikes" correspond to reads from exons, since intron and intergenic sequences generally aren't sequenced in our RNA-seq experiment. Also note how while FL and HL expression is quite high at this gene, MB looks like it has very little coverage in this region, suggesting Mecom is not highly expressed there. Scroll around to some other genes.

To avoid having to reload all the files you're visualizing each time you open and close IGV, you can save a "session", which will keep track of all the files, settings, etc. that you were using before. See more information on the IGV website: https://software.broadinstitute.org/software/igv/Sessions. 

<blockquote>
**IGV TIP**: To make things easier to visualize, you can color each track. For instance, I found it helpful to make the two replicates of each tissue type a different color. Right click on the name of the track at the left and choose "Change track color".
</blockquote>

## 5. Comparing overall expression patterns across datasets

First take a look at the `kallisto` output. You should have one directory for each experiment. For example take a look at the directory `FL_Rep1` where you should see the following files:
* `abundance.tsv`: a tab separated file giving the "TPM" values for each transcript
* `abundance.h5`: this is in a binary format (h5) that gives the bootstrap values. If you want to view this file you can use `kallisto h5dump`
* `run_info.json`: gives information about the parameters used to run kallisto

We'd like to do some sanity checks on our data. In general, we'd like to see how reproducible our results are per tissue by comparing replicates. We'd also like to do an overall comparison between tissue types. To get a rough idea of concordance between datasets, let's calculate the pairwise correlation between each one. While you could use python or R to do this, we'll use UNIX commands here for practice and since it faster to set up.

Run the following command to compute the Pearson correlation between the TPM values in both replicates of FL:
```
paste ~/week4/FL_Rep1/abundance.tsv ~/week4/FL_Rep2/abundance.tsv | cut -f 5,10 | grep -v tpm | awk '(!($1==0 && $2==0))' | datamash ppearson 1:2
```

Let's break apart this line since it introduces some new commands:
* `paste`: is a useful command to horizontally concatenate two files. Since each `abundance.tsv` file has genes in the same order, we can `paste` them together to get one big file with results from the replicates side by side.
* `cut`: is our old friend. It extracts columns 5 and 10 (which contain the two TPM columns after doing `paste`)
* We use `awk` to filter further. What does the above `awk` command do? When you present the results be sure to mention this step.
* `datamash` is used to calculate pearson correlation between columns 1 and 2. See `datamash --help` or the [datamash website](https://www.gnu.org/software/datamash/) for more info. This is a really useful tool for quickly calculating statistics about your files on the command line.

In general, to understand how a long command with many steps (such as the above) works, you can try running the command stepwise. For example first, seeing what the output of paste looks like by typing `paste ~/week4/FL_Rep1/abundance.tsv ~/week4/FL_Rep2/abundance.tsv | head` and then adding on each additional step to see what the intermediate outputs are. Try this with the command above to make sure you understand what each step is doing.

Repeat the correlation command for each pairwise analysis of all the 6 `kallisto` results. Present the results as a table or a heatmap. Which tissues were most similar? Most different? How concordant were the replicates? Are replicates more concordant with each other than with other tissues?

Make sure you save your IGV session before you leave so you can load it when we start up again on Thursday.

**That's it for today! Next time, we'll identify differentially expressed genes between the three tissues. We'll analyze regulatory elements for one of the differentially expressed genes to find specfic candidate sequences likely to be involved in limb development.**

