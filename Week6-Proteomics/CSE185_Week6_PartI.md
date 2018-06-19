# Week 6: Proteomics (part 1)
Skills covered: Mass spec analysis

This week, to give you a small break before we go into final project mode, you won't have to write a lab report. But you will have to go through the lab tutorials and fill in your lab notebook as usual. Instead of writing a report, answer the  worksheet file "CSE185_Week6_Worksheet.md" in your Github repository (bolded questions, see instructions below).

The goal of this week is to get a gentle introduction to various aspects of protein bioinformatics. Most of our work for this lab will involve exploring datasets using graphical user interfaces (GUIs). Today, you will explore some mass spec data. On Thursday, we'll get hands on experience reconstructing the sequence of a peptide from an MSMS spectra.

You are interested in studying the telomeres, which are highly repetitive regions at the ends of chromosomes. Telomeres serve many important functions, including protecting the ends of chromosomes, helping organize chromosomes in the nucleus, and allowing chromosomes to replicate properly. Each time a cell divides, telomeres shorten in length by 25-200bp. Eventually, after telomeres get too short, a cell will undergo programmed cell death (apoptosis). Thus, telomeres are an important marker for aging: generally the longer the telomeres the younger a cell.

You have performed an experiment on mouse cells to isolate all proteins associated with telomere regions. You perform tandem mass spec on the resulting mixture of proteins and would like to determine the identity of proteins in this complex mixture. 

## 1. Exploring the PRIDE database using the PRIDE inspector
Just like NCBI hosts the SRA for next generation sequencing data, EMBL-EBI (the European Molecular Biology Laboratory - European Bioinformatics Institute), hosts the PRIDE database for proteomics data from mass spectroscopy. Most proteomics data on PRIDE is stored in some version of XML format, which is not easy to read directly. So we will use PRIDE’s in-house  tool, called PRIDE Inspector, for viewing the data in a graphical user interface.

You can download the PRIDE inspector here: https://github.com/PRIDE-Toolsuite/pride-inspector.

## 2. Inspect the raw XML data file

We will be working with the file `public/week6/telomere_ms.mzXML`, which contains MS/MS spectra from proteins isolated from telomere regions. Our goal is to identify which proteins were contained in this mixture. Use:

```
less -S PATH/TO/telomere_ms.mzXML
```
to scroll a bit through this file. This is in an XML format, which can be difficult for humans to read. If you've used HTML before some of the syntax of this file might be familiar. This is a particular XML format called [mzXML](https://www.nature.com/articles/nbt1031) developed specifically for representing mass spectrometry data.

The top contains a bit of metadata about the file, and then there are individual entries for each peptide that entered the mass spec, which is reported as a scan. Because this is tandem MS data, there were two different scan levels. Level 1 reports the m/z peaks for trypsinized fragments that resulted when the peptide was first ionized to convert it to the gas phase (these fragments are often called "precursors").

In this experiment, the top 9 most abundant peaks from the first scan were subjected to additional fragmentation into smaller pieces with a process called CID (collision induced dissociation). One m/z peak at a time, each group of precursor fragments is accelerated and then the pressure is increased so that the gaseous ions collide with one another and break up. These new, smaller fragments are scanned again in Level 2, to give additional information about each original peptide.

In each reported scan, the m/z ratio and the corresponding intensity (which is proportional to abundance) of each fragment is reported as a sequence of ASCII encoded pairs.

Answer the questions below in the CSE185_Week6_Worksheet.md worksheet by manually inspecting the `telomere_ms.mzXML` file via the command line:

**What level is the first scan?**

**How many peaks are in the first scan?**

**How many scans are in the entire file?** (Hint, `cat telomere_ms.mzXML | grep "scan num="` will be helpful)

## 3. Load PRIDE Inspector and visualize spectra

Now switch over to the PRIDE inspector. Go to Open -> Open Files -> and navigate to the file `public/week6/telomere_ms.mzXML`. It may take a minute or two to load.

Click on the "spectrum" tab. The bottom window shows a list of every scan from this sample (only the first 100 are shown, you can load more in 100-scan batches with the Load Next button). The top window shows the spectra corresponding to that scan. Each peak represents an individual fragment. The y axis shows its abundance, and the x-axis shows its mass/charge ratio. Use PRIDE inspector (you may have to look around the window or try other tabs) to answer the following questions:

**In the first spectrum that loaded (spectrum ID 4970), what is the m/z ratio of the most intense peak? (You can see this by hovering your mouse over that peak).**

**In the summary tab, compute the summary charts. Then: what are the two most abundant precursor ion charges?**

## 4. Convert and filter spectra

The PRIDE inspector lets us visualize the raw data, but to understand which peptides these fragmented ions correspond to, we need to search them against a database of predicted fragmentation spectra for known proteins.

The database we will use is called mascot. It is a free, online server, but it requires a specific format and can only handle 1200 scans at time. Ideally, we would search all of our spectra, but to reduce the size of our file, we will take a subset of spectra with an intermediate number of peaks. 

We will use the tool `msconvert` from the ProteoWizard suite to convert our data to the correct format and extract only a fraction of the spectra:

```
msconvert /PATH/TO/telomere_ms.mzXML --filter "index [0,1000]" --mgf -e 300_500.mgf 
```

This filters our data to keep only scans in a specified index range. This will output a file (in a new format!) `PRIDE_Exp_Complete_Ac_31251.pride300_500.mgf`, which is a simpler format that mostly contains a list of peaks for each scan. This file can be both analyzed using PRIDE and used to search databases using the Mascot tool.

Go back to the PRIDE inspector and load our filtered `.mgf` file.

**How many spectra are in our filtered file?**

## 5. Search spetra against protein digest databses with mascot

In a web browser, go to the mascot database at http://www.matrixscience.com. Find the free search, then navigate to the MS/MS ions search. Enter:

* Your name, email, and a search title.
* Select the SwissProt database, as well as "contaminants".
* Set the taxonomy to Mus Musculus.
* Set the digestion enzyme to trypsin
* Allow up to 5 missed digestion sites
* Add the following fixed modifications: Carbamidomethyl (C)
* Add the following variable modifications: Acetyl (N-term), Oxidation (M)
* Set peptide mass tolerance to 10 ppm.
* Set Fragment mass tolerance to 0.8 Da.
* Set peptide charge to 2+ and 3+
* Check monoisotopic
* Set data format to mascot generic
* Upload your filtered `.mgf` file.
* Start search (may take 2-3 minutes)

## 6. Analyze mascot results

Your search results will include some summary information at the top, and the blue-ish chart will show you the threshold for significance based on how well your spectra agree with the listed peptide sequence.

Beneath this, in the Protein Family Summary, go to the "report builder". You will see a numbered list of protein hits, with columns giving the progein name, match score, and additional info. Click on a family number to see more details about the peptide(s) that were identified for each protein.

**The top hit is a contaminant! What is it? Give one sentence explanation of why this shows up in our list.**

**Sort the list by score. List the top 5 scoring protein hits (not contaminants!)**

**How many peptides were found to match the top hit, NUCL_MOUSE (Nucleolin)?**

**Briefly investigate the function of the top 5 genes. Are any known to be related to telomere function?** (3pts)

***That's it for today! Next time we'll learn how to reconstruct peptide sequences from tandem MS data.***


