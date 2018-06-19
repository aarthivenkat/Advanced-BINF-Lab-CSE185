import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('lab5_gwas_eyecolor.eigenvec', sep=' ', header=None)[[2,3]]
fig = plt.figure()
plt.scatter(df[[2]],df[[3]])
plt.xlabel('PCA1')
plt.ylabel('PCA2')
plt.title('Eye Color PCA Analysis')
plt.savefig('pca1_vs_pca2.png')
