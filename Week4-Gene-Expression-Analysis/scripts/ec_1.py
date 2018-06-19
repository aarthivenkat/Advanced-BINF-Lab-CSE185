import matplotlib; matplotlib.use('Agg')
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

f = open("tpm_matrix.txt", 'r')
m = []

for i,line in enumerate(f):
    col = [float(x) for x in line.strip().split('\t')]
    newcol = []
    for x in col:
        if x > 0: newcol.append(math.log10(x))
        else: newcol.append(x)
    m.append(newcol)

map = np.array(m)
plt.figure(figsize=(50,50))
g = sns.clustermap(map)
plt.suptitle("TPM_Clusters")
g.ax_heatmap.set_xticklabels(['FL_Rep1', 'FL_Rep2', 'HL_Rep1', 'HL_Rep2', 'MB_Rep1', 'MB_Rep2'])
plt.savefig('tpm_clustermap.png')
