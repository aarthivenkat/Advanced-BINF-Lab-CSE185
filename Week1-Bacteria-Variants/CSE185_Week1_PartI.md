
# Week 1: What causes antibiotic resistance? (part I)
Skills covered: basic UNIX navigation, Github, FASTQ format

Microbial resistance to antibiotics is a major problem in the treatment of bacterial infections. Evolution is a master algorithm at finding ways to survive, and there is often more than one way for bacteria to become resistant to a particular antibiotic. For scientists and doctors, knowing the exact mechanism of resistance that evolved in a bacterial strain is very important. This knowledge can guide development of new, better drugs, or even help doctors to decide which alternative antibiotics to use if one is failing to treat a patient. 

Real sequencing data from a strain of *E. coli* resistant to the antibiotic ampicillin is located in the public folder on the course server. Your job is to analyze that sequencing data to locate the mutations responsible for giving *E. coli* its antibiotic resistance property. You will research the genes that are mutated to identify the mechanism of antibiotic resistance in each case, and make recommendations for alternative antibiotics a doctor could use to treat each strain. 

But first, we have a couple of set-up items to take care of. We'll get started today by logging into the course server and getting our lab notebooks set up on Github.

<blockquote>
**NOTE: Throughout the tutorials, commands that will be entered directly into the UNIX terminal will be in a special format:**
</blockquote>

```shell
echo "Hello CSE 185"
```

## 1. Course server login and basic UNIX navigation

You can ssh in to the course server from the terminal on any internet-connected mac or Linux machine. Launch terminal, play with the preferences to chose a color scheme you like (I also like to set my cursor to blinking), then log in with the command below. Replace XXXXX with your username obtained in class. If you don't know your usename, you can look it up [here](https://sdacs.ucsd.edu/~icc/index.php). If you are not using a course-specific account, you'll need to run the `prep` command to get set up for CSE185 specific resources.

```shell
ssh XXXXX@ieng6.ucsd.edu
```

Enter the password when prompted. You will be placed into your unique home directory (`/home/[username]`). You can use
`pwd` (print working directory) to see the path to the current directory. 
```shell
pwd
```

Use `ls` (list)  to see what’s in the current directory (it should be empty). 
```shell
ls
```

Above your home directory, there is a high-level directory for the whole course. This is where all the
raw data will be located. To get to this directory, you use the `cd` (change directory) command. The
general format for this command is:
```shell
cd [directory]
```

To use the command, you replace the part in brackets (get ride of the brackets too), with the path to
the directory that you’d like to change too. (We will use a similar format throughout the tutorials for code you will need to fill in.) This path can be absolute or relative. If you just type `cd`
alone, the shell will take you to your home directory. To specify relative paths in the [directory] part of 
the command, a single period refers to the current directory and a double period refers to its parent.
So to navigate to the high-level course directory from your home directory, type:
```shell
cd ..
```

Now when you type `ls`, you should see the home directories of everyone in the class, as well as a
public directory. You only have access to your directory and the public directory, and the public
directory is read-only. Go ahead and `cd` into the public directory.

<blockquote>
**UNIX TIP**: Unix has an ‘autocomplete’ feature that will help you correctly type names and paths. If
you start typing the command below and then press the tab key, unix will automatically fill in
the rest of the directory name, and you can just hit enter. Try it. 
</blockquote>

```shell
cd public/
```

Use `ls` to see what’s in the public folder, and `pwd` to get the absolute path to the public folder. 

<blockquote>
**UNIX TIP CONTINUED**: If there are multiple options in a file that start with the same letters (ie `week1`
and `week2`), when you press tab after you start typing, the shell will autocomplete the shared part,
then beep (if the sound is on) and wait for you to specify the rest, then you can keep typing and
tabbing. 
</blockquote>

Try using autocomplete and tab to `cd` into the `public/week1` folder. Here you should see 3 files. Two of
the files are raw Illumina sequencing reads from shotgun sequencing of an E. coli strain that is
resistant to the antibiotic ampicillin. (1 and 2 refer to forward and reverse, this was a paired end run).
One file is a reference sequence of the parental (unevolved, not resistant to antibiotics) *E. coli* strain. 

While you are in the `week1` folder, compare the size of each of these files with the disk usage `ls`
command. The optional flag `-l` makes the output in "long format" and the `-h` makes the results human readable (in 
`-h` makes the results human readable (*e.g.*, 1K 234M 2G instead of number of bytes).

```shell
ls -lh
```

Answer the <font color="blue">i>Clicker question</font> about the file sizes. When you get to this point, please wait for others to catch up before starting section 2. You can start the assigned reading (see the README for this repository) if
you’d like.

<blockquote>
**IMPORTANT NOTE**: Data analysis you do in the class will all be done in your own home directory,
and in most cases, you will use an absolute path to refer to the raw data location in the public folder,
so lets just leave that data where it is, and cd back your home directory. Unless explicitly directed, **DO
NOT copy** the raw sequencing files from the public folder into your folder. This is because they are
very large and the server space for each account and the course
as a whole is limited, so we won’t make copies unless we have to. 
</blockquote>

## 2. Github setup

We will be keeping lab notebooks on github (see the lab notebooks handout for more info). Git is a framework for version control used for keeping track of changes to files and for coordinating work on files across multiple people. Github is a web-based hosting system for Git. There's a lot to learn about Git and Github. This week we'll keep it simple and learn about Git using the web interface. We'll learn more and more about Github as we go along.

To get started:

1. Make a Github account (github.com) and send your username to the TA for your section on slack.
2. From the [course homepage](https://github.com/gymreklab/cse185-spring18), click on the assignment 1 invitation link to create your own private repository for the week 1 assignment.

Your new repository, named "cse185-week1-<username>", is your own private copy of the repository for this week's assignment. Every week you'll create a similar private repository. You should see the following files listed:

* `README.md`: Every proper Github repository has a readme file with some basic info.
* `CSE185_Week1_PartI.md`: This file, which gives instructions for part I of the lab.
* `CSE185_Week1_PartII.md`: Instructions for part II of the lab.

In addition, you'll see a mostly empty file `CSE185_Week1_LabNotebook.md` where you will need to fill in your lab notebook for the week. 

Finally, you'll see a folder `labreport/` where you will upload anything used in your lab writeup for the week in addition to the final report. A template has already been started for you.

For this week, we'll work on the lab notebook and upload files using the web interface. Next week, we'll learn how to interact with github directly through the command line. We'll start by making a small edit to the lab notebook.

The tutorial and notebook files are all written using <a href="https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet">Markdown</a>, a handy language for formatting text files. You'll catch on to Markdown syntax as we go.

Start by clicking on the lab notebook file. Click the pencil icon on the top right of the file, which will allow you to edit the file on the screen. Use the editor on the screen to update the name and date at the top of the file. Click on "Preview changes" to see how the markdown will get parsed. When you're done, scroll to the bottom to commit the changes. Every time you commit, you should add a helpful message about what changes you made. Click "commit" to save the changes. 

<blockquote>
**REMINDER: Any time you enter a command that processes, measures, or analyzes data, copy that command and any relevant results into your lab notebook. See the examples in the notebook for helpful tips on how to format text and code to make your notebook readable.**
</blockquote>

## 3. Inspect raw sequencing data manually

Next you will use a software program that analyzes the overall run quality and makes a report.
Before that however, it’s good practice to manually verify a few things directly. First, you should
always open the sequence files and verify that the format is correct. You don’t want to open the entire
file, you will just do a quick inspection of the first few reads. Use the `head` command to look at the first
20 lines in your file (replace [filename.format] with the actual file name).

```shell
head -n 20 [file.format]
```

Take a look at each file. Each read has 4 lines of information, and then the next read starts on the
following line. The first line starts with the "@" symbol, and contains identifiers and information about
this read. The next line contains the actual sequence, then on line three there is a "+", which may
sometimes have the identifier and info repeated. 

Line 4 contains the quality string, where ASCII characters encode the quality score for each base.
The quality score ranges from 0 to about 40; the higher the number, the greater the confidence of the
base call. To get the actual quality score, you need to figure out the value of the symbol, then subtract
33 (this is called "Sanger" scaling type, because it's the same scaling that people use with traditional
sanger sequencing). With some older Illumina data (pre version 1.8), you subtract other numbers, like
64, for more info, see: 

http://drive5.com/usearch/manual/quality_score.html

Our data has sanger scale quality scores. Look at the beginning of the `amp_res_1.fastq` file, and find
the third read. Note its unique read identifier, and find the quality symbol for the read’s first base.
Look up the value associated with this ASCII symbol online, then answer the <font color="blue">i>Clicker question</font>.
Feel free to keep working in this section, but please listen for the correct answer.

Use `cat` to open the entire fasta reference file `NC_000913.3.fasta`. Do you notice anything different about it?

```shell
cat NC_000913.3.fasta
```

<blockquote>
**UNIX TIP**: If you accidently open a huge file and just see characters flying down your screen, matrix
style, or if you think a process you are running may be stuck, you can press control+c to cancel the
last command. 
</blockquote>

Another thing you should check directly from the shell is how many reads are in each fastq file. Some
reads will get removed during analysis, so it’s important to know what you started with. Use word
count with the lines flag to see how many lines there are in each fastq file. 

```
wc -l amp_res_1.fastq
```

From the line count, use what you know about the fastq format to calculate the number of reads in
each file, and record in your lab notebook. Additionally, answer the <font color="blue">i>Clicker question</font> about coverage before moving on.

<blockquote>
**UNIX TIP**: If you are going to be executing nearly the same command, with a small change, like you
will need to here to get the line count of the `_2` file, instead of retyping, you can press up and down at
the command prompt to scroll through your recent commands, then use the arrow keys and delete to
modify them.
</blockquote>

***That's it for this session! This was just a brief introduction to the data structure from illumina
sequencing runs. Next time, you will preprocess the data, then align our reads to the E. coli
genome and locate mutations in our antibiotic resistant strain.***

**Acknowledgements**: Adapted from a lab originally written by Dr. Katie Petrie
