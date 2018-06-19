# CSE 185 Lab 6 Worksheet

This worksheet is worth the same amount as a lab report (based on the percent of points scored for each question below scaled to a 10 point scale). 

## Part 1: Mass spectrometry

**What level is the first scan?**  
msLevel="2". Level 2  

**How many peaks are in the first scan?**  
peaksCount="205". 205  

**How many scans are in the entire file?**   
<scan num="21409". 21409 scans  

**In the first spectrum that loaded (spectrum ID 4970), what is the m/z ratio of the most intense peak? (You can see this by hovering your mouse over that peak).** (1pt)  
The m/z ratio is 555.065.  

**In the summary tab, compute the summary charts. Then: what are the two most abundant precursor ion charges?** (1pt)  
The most abundant precursor ion charges are charges 2 and 3, with frequencies 15231 and 3347, respectively.  

**How many spectra are in our filtered file?** (1pt)  
There are 1001 spectra.  

**The top hit is a contaminant! What is it? Give one sentence explanation of why this shows up in our list.** (1 pt)  
The top hit is Trypsin, which is the digestion enzyme we needed to include to perform the cleavage of cutting the C-term side of KR unless the next residue is P.  

**List the top 5 scoring protein hits (not contaminants!)** (1 pt)  
In order, the top scoring hits are NUCL_MOUSE, DHX9_MOUSE, TR150_MOUSE, EF1A1_MOUSE, PRP8_MOUSE.  

**How many peptides were found to match the top hit, NUCL_MOUSE (Nucleolin)?** (1pt)  
7 peptides, 63 amino acids (8% coverage) matched the top hit.  

**Briefly investigate the function of the top 3 genes. Are any known to be related to telomere function?** (3pts)  
1. Ncl codes for nucleolin, a nucleolar phosphoprotein involved in the production of ribosomes. According to [Khurts et al](http://www.jbc.org/content/279/49/51508.full), nucleolin interacts with telomerase and alters its subcellular localization, which in turn plays a role in telomere function. This could potentially be related to telomere function as a whole.  

2. Dhx9 codes for a helicase, implicated in transcriptional and translational regulation. [Lee et al](https://www.ncbi.nlm.nih.gov/pubmed/24990949) discussed how the suppression of the helicase may induce premature senescence, through disturbing the the replicative telomere-dependent process.  

3. Thrap3 codes for pre-mRNA splicing factors. [A recent study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4693264/#B178-biomolecules-05-02935) has revealed a relationship between DNA damage at telomeres and splicing factor expression, so perhaps there is a relation between Thrap3 and telomeres in terms of influencing stability.  

## Part 2: De novo peptide sequencing

**How many peaks are in this file? What is the range of m/z ratios? (i.e. min and max). How does the max ion mass here compare to the precursor ion mass? What is the most abundant m/z ratio?** (4pts)  
There are 1108 peaks in this file. The minimum m/z ratio is 279.199249 and the max is 1958.194306. Because the m/z ratio is equal to the mass of each ion, we know that the max ion mass here is equal to the precursor ion mass. The most abundant is 768.640625, with intensity 7950.795898.  

**Include your table in your worksheet.**(3 pts - 1 for including table, 2 for explanation of at least 2 more AAs)  

| Amino acid | Position (b-series) | Position (y-series) | Peak (b) | Peak (y) |
|--------|---------|---------|-------|------|
| I | 1 | 18 | 114.10 | 1845.93 |
| A | 2 | 17 | 185.12| 1775.89 |
| G | 3 | 16 | 305.28 | 1719.87 |
| P | 4 | 15 | 355.15 | 1623.82 |
| P | 5 | 14 | 468.37 | 1527.77 |
| S | 6 | 13 | 580.36 | 1441.74 |
| G | 7 | 12 | 677.54 | 1385.72 |
| G | 8 | 11 | 746.64 | 1329.70 |
| **P** | 9 | 10 | 822.50 | 1233.65 |
| A | 10 | 9 | 936.61 | 1163.62 |
| I | 11 | 8 | 1022.51 | 1051.54 |
| P | 12 | 7 | 1150.57 | 955.49 |
| E | 13 | 6 | 1263.66 | 826.45 |
| **Q** | 14 | 5 | 1362.66 | 699.36 |
| S | 15 | 4 | 1506.58 | 613.33 |
| M | 16 | 3 | 1603.90 | 483.29 |
| G | 17 | 2 | 1754.99 | 427.29 |
| K | 18 | 1 | 1812.90 | 145.0977 |

**Explain your answer for at least two more amino acids.**  
For the P bolded in the table, I determined this was a proline by working from the end of the b-series up, so I found a peak at 936.61 and a peak at 822.50, and I calculated 936.61 - 822.50 - 17 (C-terminus) = 97.11, which is very close to the mass of Proline, 97.05 Da. To confirm this, I found the position in the y-series by working from the top of the table down, so I found a peak at 1329.70 and calculated the position of the next peak if this amino acid were truly proline as 1329.70 - 97.05 + 1 (N-terminus) = 1233.65. It turns out that there is a true peak close to 1232, so this reified my answer.  

For the Q bolded in the table, I went through a similar process of working backwards to find the amino acid from the b-series, and then checking my answer by seeing there exists peaks where I anticipate in the y-series. 1506.58 (peak 1) - 1362.66 (peak 2) - 17 = 126.92, which is close to 128.06 Da, the mass shift for Q, glutamine. 826.45 - 128.06 + 1 = 699.36, the putative position in the y-series, and it turns out that there is a true peak near 695, so I can confirm this amino acid.  


**This was pretty tedious! If you were going to write a program to automate this process, how would you do it? (max 2-3 sentences about your idea.** (2 pts)  

If I were to code this, I would set a threshold for peaks, smooth the peaks so that there are obvious intervals between peaks, and once I find the difference between the two peaks in both series, determine if the value is close enough to a mass shift for an amino acid to implicate the amino acid at that position. Further, I would attempt to use optimization and high resolution datasets to be able to accurate assess the amino acid when the mass shifts are very close (e.g. isoleucine and leucine).  

**What is the op scoring sequence? How close is it to your guess? Note: neither answer is probably completely correct!** (1 pt)  
The top scoring sequence is PQELLSQLVQ. My guess was IAGPPSGGPAIPEQSMGK. This is not close whatsoever.  

**Hypothesize why some peaks are missing from our spectrum file.** (1 pt)  

Correct:  
![](https://github.com/cse185-sp18/cse185-week6-aarthivenkat/blob/master/FragIon.PNG)  

I can infer that the procedure against which the mass spec experiment requires various calibrations to optimize the data, such as during the separation before MS and with the MS instrument. Further, perhaps the MS instrument itself was at low resolution, or, because peaks are defined by their relative intensity to the base (largest) peak, some peaks may have had a high intensity but a low relative intensity.  

**Suggest a way we could modify our mass spec experiment to distinguish between amino acids with near identical masses, such as Isoleucine and Leucine.** (1 pt)  
We could use a higher resolution MS instrument, which would be more expensive but could resolve masses to more decimal places.  
