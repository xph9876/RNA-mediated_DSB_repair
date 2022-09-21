import sys
import os
sys.path.append("C:\\Users\\Tejasvi\\Code\\DSB_Project\\RNA-mediated_DSB_repair\\analyze_nhej\\utils")
import pandas as pd
import glob

freq = []
for x in glob.glob("libraries_3/*"):
  print(x)
  if ('30bpDown' not in x) and ('noDSB' not in x):
    dat = pd.read_csv(x, sep="\t")  
    for col in dat.columns[dat.columns.str.startswith("Freq")]:
      freq.append(dat.loc[dat[col] <= 1e-5, col].sum())

print(max(freq))

x = pd.read_csv("analyze_nhej/layouts/universal/2DSBanti_CD/sequence_data_withoutSubst.tsv", sep = "\t")