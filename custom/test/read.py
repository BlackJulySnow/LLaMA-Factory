import os
import pandas as pd
from rdkit.Chem import BRICS
from rdkit import Chem
from tqdm import tqdm
base = "C:\\Users\\B1GGersnow\\Desktop\\mol"
def brics(s):
    mol = Chem.MolFromSmiles(s)
    return ".".join(BRICS.BRICSDecompose(mol))
smiles = []
for dir in os.listdir(base):
    d = os.path.join(base,dir)
    if not os.path.isdir(d):
        continue
    for f in os.listdir(d):
        df = pd.read_csv(os.path.join(d, f), sep=' ')
        smiles.append(df)
result = pd.concat(smiles, ignore_index=True)
result = result.drop_duplicates(subset='smiles')
tqdm.pandas(desc="Converting SMILES to BRICS")
result['brics'] = result['smiles'].progress_apply(brics)

result.to_csv('output.csv', index=False)