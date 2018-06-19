# Week 6: Proteomics (part 2)
Skills covered: More mass spec, protein sequencing

On Tuesday, we were introduced to mass spec data and different tools for manipulating it. Today, we'll dig a bit deeper and learn how to reconstruct the sequence of a short peptide given it's spectrum.

You have performed one round of mass spec and identified the most abundant component (the precursor ion) which you would like to study further. You subject the ion to CID in order to fragment the peptide and perform an additional (tandem) round of mass spec to obtain its spectrum.

## 8. Examining the spectrum

The spectra for this peptide is given in `public/week6/example_peptide.mgf`. We were introduced to `.mgf` format on Tuesday. Use `head -n 20` to take a look at this file. You should see:

* The top says `BEGIN IONS`. The last line says `END IONS`. These act as bookmarks to tell where one spectrum ends and the next begins in the file. Here we only have one spectrum so you only see these once.
* The next line tells us the mass of the precursor ion, which in this case is 1958.
* The third line gives the charge. The charge here is just +1, which is convenient since now all of our m/z ratios are just equal to the mass of each ion.
* The next lines give the peaks of the spectrum. Column 1 gives the m/z ratio (which again, is the same as mass here since z=1). Column 2 gives the intensity.

**How many peaks are in this file? What is the range of m/z ratios? (i.e. min and max). How does the max ion mass here compare to the precursor ion mass? What is the most abundant m/z ratio?**

Now, load this file into the PRIDE inspector tool we were using on Tuesday. Refer back to Part I instructions for how to run this tool.

Go to the "Spectrum" tab. There you should see the peaks represented in the `.mgf` file. Do your answers to the previous question agree with the visualization?

## 9. Reconstruct the pepdtide sequence

We will use the spectrum to reconstruct the sequence of the precursor ion. As a hint, our peptide is 18 amino acids long. As you figure out different amino acids in the forward and reverse direction it is a good idea to keep track of them in a table somewhere.

As a refresher, refer to the diagram below:

![peptide sequencing](https://github.com/gymreklab/cse185-spring18-week6/blob/master/peptide_seq_example.png)

We will walk through the spectrum (from right to left, or heaviest to lightest mass) to determine which amino acid is represented by the distance between each pair of peaks. Recall that after CID, our peptide gets broken up into both "y" type ions (which extend from the C-terminus) and "b" type ions (which extend from the N-terminus). So our resulting spectrum will have peaks corresponding to the mixture of these two series which we'll have to deconvolute: the y-series has peaks for each amino acid read in one direction, and the b-series has peaks reading in the opposite direction. 

Before we get started, notice:
* There are a ton of peaks of varying abundance. Many of these are noise! In general, the biggest spikes (highest abundance) are real peaks and the tiny ones might not actually represent real pieces of our peptide.
* Some peaks will be missing altogether. So as we put together this puzzle, you'll have some question marks we'll have to fill in later.

We started with a peptide (the precursor ion) with mass 1958. In PRIDE, zoom in on this region to get a better look at the peaks in this region (which are much smaller and thus difficult to see when we are completely zoomed out). In general the peaks for the highest mass peptides are pretty low abundance. It might be easier to look at the `mgf` file to better see where some of the lower abundance peaks are.

To reconstruct the peptide, let's start at the "parent mass" of 1958 and work back from here. We'll use the following table of masses for each amino acid:

| Mass shift (Da) | Name | Abbreviation |
|-----------------|-------|------------|
|57.021464|Glycine|G|
|71.037114|Alanine|A|
|87.032028|Serine|S|
|97.052764|Proline|P|
|99.068414|Valine|V|
|101.047679|Threonine|T|
|103.009185|Cysteine|C|
|113.084064|Isoleucine|I|
|113.0841|Leucine|L|
|114.042927|Asparagine|N|
|115.02694|Aspartic acid|D|
|128.058578|Glutamine|Q|
|128.094963|Lysine|K|
|129.04259|Glutamic acid|E|
|131.040485|Methionine|M|
|137.05891|Histidine|H|
|147.068414|Phenylalanine|F|
|156.10111|Arginine|R|
|163.063329|Tyrosine|Y|
|186.079313|Tryptophan|W|

Note some have very similar masses! e.g., Isoleucine and Leucine will be pretty impossible for us to tell apart! So we won't be able to completely resolve our sequence.

In addition to these masses, we'll have to adjust for a mass of 1.0078 for the N-terminal amino acid (farthest right amino acid in the y-series, farthest left in the b-series) and 17.0027 for the C-terminal amino acid (fatherest left in the y-series, farthest right in the b-series).

So, starting with 1958, we can calculate all possible shifts for the y and b series based on what the first (N terminal) or last (C terminal) amino acid is in our sequence. Let's do the first one as an example. To get the list of possible shifts for the y-series, we'll take 1958-mass+1.0078 for each possible amino acid. For the b series, we'll take 1958-mass-17.0027 for each. You should get:

| Amino acid | Expected peak (y series) | Expected peak (b series) |
|-------|-----|------|
|G|1901.986336|1883.975836|
|A|1887.970686|1869.960186|
|S|1871.975772|1853.965272|
|P|1861.955036|1843.944536|
|V|1859.939386|1841.928886|
|T|1857.960121|1839.949621|
|C|1855.998615|1837.988115|
|I|**1845.923736**|1827.913236|
|L|1845.9237|1827.9132|
|N|1844.964873|1826.954373|
|D|1843.98086|1825.97036|
|Q|1830.949222|1812.938722|
|K|1830.912837|**1812.902337**|
|E|1829.96521|1811.95471|
|M|1827.967315|1809.956815|
|H|1821.94889|1803.93839|
|F|1811.939386|1793.928886|
|R|1802.90669|1784.89619|
|Y|1795.944471|1777.933971|
|W|1772.928487|1754.917987|

Take a look at the `.mgf` file to see if you have peaks corresponding to any of these masses. Looking at the tail (e.g. `tail -n 30 example_peptide.mgf`) you'll notice a peak around ~1813. That is close to the expected peak for lysine (K) in the b-series! So we can guess that the last amino acid in our peptide is K (which makes sense, this peptide resulted from trypsinization, which results in peptides ending in "K" or "R"). For the y-series, we actually should have gotten a peak at ~1845, since the first amino acid in our peptide is I (you already have 2/18 figured out!). That peak isn't in our data though... (biology is messy!).

So we can already start to fill out our table:

| Amino acid | Position (b-series) | Position (y-series) | Peak (b) | Peak (y) |
|--------|---------|---------|-------|------|
| I | 1 | 18 | | 1845.923 |
| ? | 2 | 17 | | |
| ? | 3 | 16 | | |
| ? | 4 | 15 | | |
| ? | 5 | 14 | | |
| ? | 6 | 13 | | |
| ? | 7 | 12 | | |
| ? | 8 | 11 | | |
| ? | 9 | 10 | | |
| ? | 10 | 9 | | |
| ? | 11 | 8 | | |
| ? | 12 | 7  | | |
| ? | 13 | 6 | | |
| ? | 14 | 5 | | |
| ? | 15 | 4 | | |
| ? | 16 | 3 | | |
| ? | 17 | 2 | | |
| K | 18 | 1 | 1812.90 | |

We can also go ahead and calculate where we'd expect peaks starting from the other direction. If the first amino acid is "I", for the b-series we'd expect a peak at that mass (113.0841+1.0078 since it's the N-temrinus). If the last amino acid is "K", for the y-series we'd expect a peak at the mass (128.0950+17.0027). (We don't have m/z ratios in that range but keeping track of these peak values can be helpful for finding the next offsets).

| Amino acid | Position (b-series) | Position (y-series) | Peak (b) | Peak (y) |
|--------|---------|---------|-------|------|
| I | 1 | 18 | 114.0919 | 1845.923 |
| ? | 2 | 17 | | |
| ? | 3 | 16 | | |
| ? | 4 | 15 | | |
| ? | 5 | 14 | | |
| ? | 6 | 13 | | |
| ? | 7 | 12 | | |
| ? | 8 | 11 | | |
| ? | 9 | 10 | | |
| ? | 10 | 9 | | |
| ? | 11 | 8 | | |
| ? | 12 | 7 | | |
| ? | 13 | 6 | | |
| ? | 14 | 5 | | |
| ? | 15 | 4 | | |
| ? | 16 | 3 | | |
| ? | 17 | 2 | | |
| K | 18 | 1 | 1812.90 | 145.0977 |

Your goal in this section is to manually reconstruct as much of the sequence as you can. You can approach this in any way you want. A couple suggestions:

* You can walk through the peaks in either direction. For example, you can look at the `.mgf` file and find the next peaks. Starting from the tail, there seem to be peaks at ~1776 and ~1755. Those are probably the next peaks in the y- and b-series. Figure out what amino acids those might correspond to, and keep walking backwards.
* You could also try to first identify the major peaks, then calculate all shifts of nearby peaks to figure out what amino acids those shifts correspond to.

Try to construct as much as you can. You may work with a partner. One of you might want to work from left to right and the other from right to left to go through the peaks and identify likely next amino acids. Fill in the table above as you go. 
Fill in as much as you can, focusing on understanding the method. Spend at most ~45 minutes on this so you have time for the rest. **Hint, if you get stuck, the true sequence is given in part 11 below. But try to reconstruct some of the sequence yourself first.**

**Include your table in your worksheet (with at least 2 more amino acids filled in). Explain your answer for at least two more amino acids.**
 
**This was pretty tedious! If you were going to write a program to automate this process, how would you do it? (max 2-3 sentences about your idea.**

## 10. Programmatically determine the peptide sequence

Now, we'll use a program called `pepnovo` (developed here at UCSD!) to programmatically determine the peptide sequence. This tool uses a probabilisitc network to model peptide fragmentation events from mass spec. Run `pepnovo` on our spectrum with the following command:

```
pepnovo -file public/week6/example_peptide.mgf -model_dir public/week6/Models -model CID_IT_TRYP 
```
where you will of course have to slightly modify the paths to the `public` directory. This command took in:

* `-file` gives the `.mgf` file containing all of our peaks
* `-model CID_IT_TRYP` tells it to model fragmentation based on performing CID on trypsinized fragments, which is how our spectrum was obtained.
* `-model_dir` tells it to look in a non-standard location for the model files

This will output candidate peptide sequences ranked by score.

**What is the top scoring sequence? How close is it to your guess? Note: neither answer is probably completely correct!**

## 11. Compare to the optimal spectrum.

The original sequence of the peptide is actually IAGIIKVIKSLIEQFTGK. Go to this site: http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html to determine the theoretical spectrum for this peptide. Enter the peptide sequence, choose "mono", "+1", "b", and "y". Compare the results to your table above.

**Hypothesize why some peaks are missing from our spectrum file.**

**Suggest a way we could modify our mass spec experiment to distinguish between amino acids with near identical masses, such as Isoleucine and Leucine.**

## 11. Finish the worksheet

When you've completed the worksheet, make sure your repository is up to date on Github (by doing `git push`). This is due at 11am next Tuesday.

***Acknowledgement: Jeramie Watrous conceived of this assignment and provided MS data.***
