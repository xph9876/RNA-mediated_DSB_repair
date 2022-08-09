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