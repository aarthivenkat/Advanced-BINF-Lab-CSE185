
# Week 1: What causes antibiotic resistance? (part II)
Skills covered: fastqc, data filtering, sequencing alignment, variant calling, variant effect prediction

Today we will start by taking a look at the quality of our sequencing data, and carry out some
preprocessing steps. 

## 4. Inspect raw sequencing data with fastqc

Cd into your home directory (type `cd`). Make a directory that will contain all of the data for week1, use
`ls` to confirm its there, and `cd` into it. 

```shell
mkdir week1
```

This is where you will be working for the rest of the week. 

OK, now you are ready to run `fastqc`, a simple fastq statistics analysis program. First, make sure the
program is working and properly installed by typing the command below. You should see the manual
page. Please raise your hand if this is not working.

```shell
fastqc -h
```

Run the program fastqc on the two fastq files. You will have to specify the path to the `public/week1`
folder, since you should leave the files there. You also have to tell fastqc to `-o` output the files to the
current directory, which you indicate with a `.` The command to do all this is below, but you have to
specify the full root path to each fastq file (*i.e.* `/home/linux/ieng6/cs185s/public/week1` )

```shell
fastqc -o . [root path]/amp_res_1.fastq [root path]/amp_res_2.fastq
```

Check with `ls` that this generated some files. The html files contain the whole report. To look at them,
you will need to download them to your local computer. 

Open a new terminal window and use `pwd` and `ls` to find the path to the desktop. Then use "secure
copy" with the command below, using the full rooted path to
the html file and your local desktop. If you’re not sure of the path, just `pwd` while you're in
the relevant directory. Copy both files. 

```shell
scp [username]@ieng6.ucsd.edu:[root path]/file.html [path to desktop]
```

You will be prompted for your cs185 password. Once the files appear on your desktop, open them
and take a look. Do the basic statistics match what you calculated for the number of reads last
time? On the left, you’ll see a navigation window with green (normal), yellow (slightly abnormal), and
red (very unusual) circles for several kinds of data analysis. If you have any red circles, record them
in your notebook, and read the fastqc documentation on the analysis modules to try to learn what
they mean. 

http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

**For you lab report**, upload these files to your `labreport` folder on github so you’ll have them when you write your lab report.
What do you think we should do about anything fastqc identified as unusual? Answer the <font color="blue">i>Clicker</font> question about it, then feel free to move on. 

## 5. Filtering the reads

Removing low quality data can improve our downstream analyses. Look at the results of your `fastqc` analysis under the category "Per base sequence quality." Notice that the quality scores tend to decrease toward the end of the reads. Trimming these low quality bases is one way to remove sequencing errors.

There are a large variety of programs designed to improve the overall quality of your sequencing
reads before proceeding with downstream analysis. We are going to use one trimming program, called sickle.

Sickle uses a sliding window trimming algorithm, where a small window slides along the sequence, starting from one end, cutting off sequence until the average quality is above a user specified threshold, then the window keeps going until the average quality drops below that threshold, and cuts off any sequence past that point. 

To learn
about this command, just type it as shown below, without any input files. 

```shell
sickle pe
```

You can specify both the quality threshold and a length threshold, which will discard any sequences that are shorter than a certain length after trimming. The defaults for both of these values are 20. Run sickle on your data with the default settings, remembering to specify the full paths to the fastq files in the public folder. Below is an example command.

```shell
sickle pe \
    -f [root path]/file1.fastq \
    -r [root path]/file2.fastq \
    -t sanger \
    -o trimpair1.fastq \
    -p trimpair2.fastq \
    -s singletons.fastq
```

Note this command should all be on one line. To make it more readable it is broken up into multiple lines. A `\` means the command is continued on the next line. Remember, 1 and 2 here are our forward and reverse reads, located in the public file. The general format is the command (`sickle pe`) followed by flags that tell sickle what the file that comes after the flag is. 

With paired-end data alignment algorithms often use the fact that matching forward and reverse reads
from the same DNA molecule can’t be that far apart (because there is an ideal molecule size for
cluster generation in the flowcell). Once a forward read is mapped, the aligner knows it doesn’t have
to scan the reference sequence very far to find where the corresponding reverse read (from the same
molecule) could go. Because aligners use this information, it’s very important to make sure that the
number of reads in the forward and reverse files stay coordinated. There are MANY trimming
programs out there for next generation sequencing data, but only some that can handle paired end
data properly. If you use something other than sickle in the future, be sure you choose carefully and always check to make sure the output makes sense! 

You’ll notice that there are 3 output files. The first two are the trimmed forward and reverse reads,
where BOTH the forward and its matching reverse read passed the trimming filter. The singletons file
contains reads where one read passed the trimming filter, but its partner did not. The reads that do
not pass are not output.

Sickle will report back stats, on how many paired reads were kept, but you should check the line
count (`wc -l`, dived by 4) of the trimmed paired files manually too. Record this in your notebook. 
To see how trimming affected the overall quality of your data, repeat the fastqc analysis you did in
section 3, but this time on the `trimpair#.fastq` files. (Note, since these are now in your working
directory, you won’t have to type out the full path when you run fastqc.)

Like before, transfer the html files for each set (1 or 2) to the desktop for viewing (this will still need
the full path), and **for your lab report, upload them to your `labreport` directory on Github.**

What happens if we increase the default quality score to 30? Try it with the command below (be sure
to name them something distinct like trimpair30_1.fastq). 

```shell
sickle pe \
    -q 30 \
    -f [root path]/file1.fastq \
    -r [root path]/file2.fastq \
    -t sanger \
    -o trimpair1_30.fastq \
    -p trimpair2_30.fastq \
    -s singletons_30.fastq
```

How many paired reads did sickle keep this time? You don’t have to redo the fastq analysis just
answer the <font color="blue">i>Clicker question</font>. 

## 6. Aligning sequences to reference

To make sense out of these small sequences of DNA, we will map them to our reference sequence
(the already-solved genome of normal, non-resistant *E. coli*). Mapping works by taking each read and
trying to align it to every possible 100 bp window in the reference. The read is mapped to the position
with the best possible alignment, indicating that the read probably corresponds to that part of the *E.
coli* genome. 

There are many alignment programs ("aligners") available. The earliest alignment algorithms (Smith Waterman
and Needleman-Wunsch) are still used today to compare small pieces of DNA one-by-one,
but the computing power needed to map millions of short reads to genomes that are millions to
billions of basepairs long requires special algorithms to speed up the process. 

Today we will use an aligner called BWA-MEM, which makes use of the Burrows-Wheeler transform
for reversible data compression (the reference is summarized with a special reversible index, the
index makes its faster to search). BWA-MEM is optimized for ‘long’ next-generation sequencing reads
of 100bp or more, which may contain several mutations, insertions, or deletions, as it works by finding
the ‘maximum exact match’ within a read to the reference, rather than forcing the entire read to
match. 

Selecting the right aligner for your data can be challenging, and there are always new ones coming
out, so it’s important to keep up with the latest developments and pay attention to benchmarking
studies (where researchers compare the performance of different algorithms on the same data set).

First, you need a copy of the reference sequence in your working directory. This is one of the only
situations where you should copy the data from the public folder for the class. First, cd into your
week1 working directory, then use the copy command below, but replace the path with the correct
path to the week1 reference in the public folder, and the path to target directory with a single "." for
your current directory.

```shell
cp [root path]/NC_000913.3.fasta [target path]
```

Next, we have to index the reference file. There is a command in bwa to do this. To see the
commands available in bwa, just type `bwa`. To see the options and usage for a specific command,
type `bwa` and that command. Do this to see how bwa index works.

```shell
bwa index
```

**Based on the usage details you see**, run bwa index on the reference sequence with the default
options. Record the command you used in your lab notebook. Check that the command worked
with `ls`. You should see a bunch of new files based on the reference.

The alignment command takes inputs in this order: first, specify the actual reference name in fasta
format (just the original reference sequence, as long as you’re in the same directory that you indexed
in, bwa will find the new index files). Then specify your actual data in fastq format. (We have paired
end data so there are two fastq files, separated by a space). Align your trimmed, paired sequences
reference with the command below (all one line). As always, replace the sample filenames below with
your ACTUAL file names. **Record the command you used in your lab notebook.**

```shell
bwa mem NC_000913.3.fasta trimpair1.fastq trimpair2.fastq > trimpaired.sam
```

<blockquote>
 **UNIX TIP**: The ‘>’ symbol is called a redirect. If we didn’t have it here, bwa mem would just output
the results to the standard out (aka the screen) which would take up way too much space AND we
wouldn’t have a record. Redirection takes the output of the preceding command and places it into a
new file specified after the ‘>’. 
</blockquote>

**Alignment may take a few minutes, in the mean time, read about "sam"**: BWA outputs data in the "SAM" format. We will dive right in with a manual inspection of the data, but to learn more about sam,
check out: https://samtools.github.io/hts-specs/SAMv1.pdf

## 7. Meet SAM, inspect the alignment

Look at first five lines of your new sam file with head:

```shell
head -n 5 trimpaired.sam
```

How many separate lines are in the output from head? How many reads? Answer the <font color="blue">i>Clicker
question. </font>

Not all reads will be successfully aligned to the reference. The sam file contains all reads, whether
they successfully aligned or not. For bioinformatics pipelines, it's important to know what fraction of
your reads aligned. If there are a lot reads that failed to align, that could indicate that your DNA was
contaminated with some other source, or that something went wrong with the alignment. 

In the SAM format, each line (after the headers) starts with a unique read ID (name.number) in the
first "field". Fields are delimited (separated) regions of data, similar to columns on a chart. Sam uses
tabs as the delimiters. In later classes we might get into what all of the different fields mean (you
should be able to recognize the quality string and the base string).

For now we are going to use the built in utilities of samtools, a program designed to read and parse
sam files, to decode the samfile for us. Run the commands below to get some basis statistics. 

```shell
samtools flagstat trimpaired.sam
```

Does the "in total" match your read count after trimming? What percentage are mapped?
These are the two most important things here. **Record these in your lab notebook. **

Next, you need to compress and sort the sam file with the commands below. A compressed sam file
is called a bam file. Like the reference, bam files need to be indexed. This gets the data ready for
some of the next commands we will use. 

```shell
samtools view -S -b trimpaired.sam > trimpaired.bam
```

```shell
samtools sort trimpaired.bam > sortedtrimpaired.bam
```

```shell
samtools index sortedtrimpaired.bam
```

Next you will visualize what the data actually looks like by using a tool called `samtools tview`, that
converts the complicated fields in the sam format into a conventional ‘alignment’ where each read is
mapped to its position on the reference. 

```shell
samtools tview sortedtrimpaired.bam NC_000913.3.fasta
```

This will bring up a simple representation of the alignment. Press the ? key to see your options for
navigation, or use the arrow keys. The reference sequence will be at the top, along with numbers
indicating the genome position, and individual reads will be shown below, aligned with the part of the
reference they map to. When the reads from your ampicillin resistant bacterial strain have the same
base as the reference there is a dot or a comma, and when they differ, you will see the base.
Press escape to exit and enter the command again. How do you interpret the A at position 46? Answer the <font color="blue">i>Clicker question.</font> To navigate to position 46 (or any other position) quickly, press the `g`
key to open Goto dialog box. Enter a position with `NC_000913.3:Position` format, where `NC_000913.3` is the name of the genome. This name can be extracted from the first line of the reference fasta file with the following command.

```shell
head -n 1 NC_000913.3.fasta
```

## 8. Make a pileup, call variants

The goal now is to go through our data, and for each position in the reference genome, see how
many reads have a mutation at the same position. The tview image looks great, but would take too
long to go through all 4 million bp of the *E. coli* genome. The SAM file is inconvenient for this, because
the reads aren’t matched with the reference in a simple way (though they do contain all the
information needed to make the alignment in tview). The solution is to make an intermediate file type
called an mpileup, because it goes through each position and “piles up” the reads, tabulating the
number of bases that match or don’t match the reference. 

Like tview, mpileup requires a sorted, indexed bam file. Run the basic command below. It may take a
few minutes.

```shell
samtools mpileup -f NC_000913.3.fasta sortedtrimpaired.bam > my.mpileup
```

See what the pileup looks like with ```head -n 100```. This file is much easier to read, but it still contains a lot
of boring information (*i.e.*, the first 100 reads mostly don’t have any mutations). We could write a script
to only spit out the positions where there are non reference bases, but one general guideline to follow
in bioinformatics is don’t reinvent the wheel if you don’t have to! There will be plenty of times when
you need to develop new parsing code or new algorithms, so in the mean time, if someone has
already written a program to do what you want, use it! 

In this case, we will be using a program called VarScan (variant scanner). Until now, all the
programs you have been using (BWA, fastqc, samtools...) were pre-installed for you in our public
folder. Here you’ll get to try your hand at installing one type of program file. 

Use the command below (curl - it's a linux/unix command to access websites) to download VarScan
into your week1 folder. Notice that we are redirecting the contents of the download into a .jar file. 

```shell
curl -L https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download > VarScan.jar
```

Try to bring up the man page with:

```shell
java -jar VarScan.jar
```

This install method won’t be the same for all bioinformatics programs, and it might not even be the
same for other java programs, as installing open source software often has dependencies
(prerequisite software) and can be sensitive to where its installed and how it is executed. We'll try to touch on this throughout the course.

On the VarScan man page, you should see several commands. We are interested in the `mpileup2snp`
command. Go ahead and enter it with `-h` to bring up the manual. 

```shell
java -jar VarScan.jar mpileup2snp -h
```

There are lots of ways we could filter our data. VarScan lets the user define their own cutoffs for
including data in calling variants, other programs rely on complex statistical procedures to evaluate
the likely hood of real mutations.

The first option we are interested in changing today is the `--min-var-frequency` option. This sets the
minimum % of non-reference bases at a position required to call it a mutation in the sample. Answer
the <font color="blue">i>Clicker question</font> about it. 

Try running the program with the threshold we decide is best (replace the N with the decimal value
corresponding to the percent we chose, ie 50% = 0.50). The `--variants` flag tells VarScan to only
output positions that are above our threshold, and `--output-vcf 1` option tells it we want the output in
yet another kind of data format called vcf (variant call format). The command may take a few minutes. 

```shell
java -jar VarScan.jar \
    mpileup2snp my.mpileup \
    --min-var-freq N \
    --variants \
    --output-vcf 1 > VarScan.vcf
```

Go ahead and take a peek at the file with:
```shell
cat VarScan.vcf
```

You should see a list of variants with a bunch of header info explaining the different fields. The key
values for us are the position, the reference, and the alternative base. These are the mutations in our
antibiotic resistant *E. coli* strain! **Record them in your notebook**, along with their "variant allele
frequency" (use the header to try to find where that is, and ask for help if you can’t find it). 

## 9. Variant effect prediction

The next task (and the final one for lab), is to find out where these mutations are, and whether they
actually change any proteins in the host (mutations can also occur outside of genes in non-coding
regions, or they can be synonymous, where they substitute a codon for they same amino acid).

You COULD go to NCBI, look up each of the variant positions in the reference sequence, see if there
is a gene there, download the gene sequence, find your variant (again), then translate the sequence
to see if the variant changes the protein sequence.

Instead, we are going to use a convenient web app that will do this for us, which can be found at
http://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_mg1655/Tools/VEP.

Unfortunately, the Ensemble Variant Effect Predictor doesn’t like the way the chromosome field is
labeled in our vcf file, so we have to edit it first to match.

For that, we will use awk. Awk is a built-in linux/unix tool. Awk excels at parsing delimited data, where
you have a lot of fields, like in a vcf file.

The simplest form of an awk command is (all one line):

```shell
awk '/search_pattern/ {actiontotakeonmatches; another_action;}' [path to file]
```

Awk has a lot of cool features, and its worth learning more. For now, we will learn by example as we
need particular awk one-liners. You will be running the command below:

```shell
awk '{if (NR>24) $1="Chromosome"; print}' VarScan.vcf > mymodVarScan.vcf
```

Here, we aren’t using a search pattern because we want to keep ALL of the lines, but INSIDE the
brackets, we’ll use an if statement to limit our action to only the lines after the header (NR stands for
number of read lines, so its saying ‘only if the line number is >24, do what’s next’. Awk uses the dollar
sign to mark fields, so what’s next is ‘set field 1 to read “Chromosome” instead of what’s there. 

Go ahead and run it, then check that the new file looks ok. 

Now that we have the correct file, use scp from a new, non ssh terminal window to put it on your
desktop. (all one line).

```shell
scp <username>@ieng6.ucsd.edu:[rootpath]/mymodVarScan.vcf [path to desktop]
```

Go to the Variant Effect Predictor website, and click on ‘new job’. Upload your file, then check the
filter that says “Return results for variants in coding regions only”. Run the job (it should take
only a few minutes), then take a look at the results! 

**Record** whether each mutation occurs in a gene, whether it is missense (changes the amino
acid sequence), non-sense (introduces a frameshift or early stop codon), or synonymous (no
amino acid change) and record what that **gene name** is. If there is an amino acid
substitution, **record** what it is. 

## 10. Write your lab report

Hoorayy!!!! You’ve made it to the end. But you’re not done yet! Simply knowing the identity of
these mutations won’t help the doctor decide how to treat a patient.

For every mutation you found that changes the protein sequence, you must research each
gene by name to find out what it does (Ecoli Wiki, EcoCyc, google, and pubmed are good
resources). Try to track down how that function could be involved in antibiotic resistance. If
one of the proteins is a gene regulator (transcription factor), try to find out what it regulates.
Include what you find in your lab report, and cite where you found it. Include this in the results
section.

Also in your results section, make a table showing how many reads you started with, how
many were left after trimming, and how many aligned. Please also include a description or
image of the per-position read quality before and after trimming.

It would require additional biochemical testing to verify EXACTLY how each mutation changes
the function of the protein its in. However, your job is to try to come up with a working
hypothesis that could explain the mechanism of resistance behind the mutations (for instance,
if the mutation is in a gene that makes a protein that is the antibiotic’s target, you could
propose that that mutation changed the target so the antibiotic couldn’t bind anymore). Use
what you know about the 4 mechanisms of antibiotic resistance to make predictions. Not all of 
the proteins may be involved in our E coli strain’s resistance, so you only need to make
predictions for THREE of the mutations, but you should research all of them, so you are using
the ones you can best make a logical prediction for. Include your predictions for the
mechanism of antibiotic resistance for three of the genes in your discussion, and explain the
logic behind your predictions.

Finally, make a treatment recommendation for someone infected with this strain of E coli. Suggest alternative antibiotics with different targets, and/or secondary therapies that might be
useful.

What do you think you learned from whole genome sequencing that traditional antibiotic
resistance testing wouldn’t have told you?

Please follow these specific Week 1 directions **in addition to** the general format guideline for
all lab reports.

**Acknowledgements**: Adapted from a lab originally written by Dr. Katie Petrie
