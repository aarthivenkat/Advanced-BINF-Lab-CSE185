#!/usr/bin/env python
"""
Script for CSE185 Lab 5 to plot GWAS results
Modify this script to make QQ plots and Manhattan plots for your GWAS results!

Usage:
./gwas_plotter.py <plink.assoc file> <outprefix>

Plink association file is a fixed width file as output from --logistic, --linear, or --assoc
Output files:
  <outprefix>_manhattan.png
  <outprefix>_qq.png
"""

import sys
try:
    datafile = sys.argv[1]
    prefix = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Set up for plotting without display enabled 
import matplotlib as mpl
mpl.use('Agg')

# Import other libraries
from matplotlib.colors import hex2color
import numpy as np
import pandas as pd

############# Manhattan plot #############
# Set up figure sizes in matplotlib
from assocplots.manhattan import *
mpl.rcParams['figure.dpi']=150
mpl.rcParams['savefig.dpi']=150
mpl.rcParams['figure.figsize']=7.375, 3.375

# Set up chromosome names and colors
chrs = [str(i) for i in range(1,23)]
chrs_names = np.array([str(i) for i in range(1,23)])
chrs_names[1::2] = ''
colors = ['#1b9e77', "#d95f02", '#7570b3', '#e7298a']
# Converting from HEX into RGB
colors = [hex2color(colors[i]) for i in range(len(colors))]

# Load data
data = pd.read_fwf(datafile)

manhattan(data['P'], data["BP"], data["CHR"].apply(str), '', \
          type='single', \
          chrs_plot=[str(i) for i in range(1,23)], \
          chrs_names=chrs_names, \
          cut = 0, \
          title='Eye color', \
          xlabel='chromosome', \
          ylabel='-log10(p-value)', \
          lines= [], \
          colors = colors, \
          scaling = '-log10')
plt.savefig('%s_manhattan.png'%prefix, dpi=300)

############# QQ plot #############

from assocplots.qqplot import *
mpl.rcParams['figure.dpi']=100
mpl.rcParams['savefig.dpi']=100
mpl.rcParams['figure.figsize']=5.375, 5.375

qqplot([data["P"]], 
       ['eyecolor'], 
       color=['b'], 
       fill_dens=[0.2], 
       error_type='theoretical', 
       distribution='beta',
       title='')

plt.savefig('%s_qq.png'%prefix, dpi=300)

