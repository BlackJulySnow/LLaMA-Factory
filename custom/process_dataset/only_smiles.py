import pandas as pd

df = pd.read_csv("custom/smiles/smiles.csv", usecols=["Smiles"], low_memory=False)
df.to_csv("custom/smiles/save.csv",index=False, encoding="utf-8")