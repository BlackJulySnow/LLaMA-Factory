import os
import json
from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import MACCSkeys
from tqdm import tqdm

with open("custom\dataset\ChemDual_fragment_test_brics1.json") as f:
    data = json.load(f)

step = [1,3,5,10]
fig1 = []
fig2 = []
result1 = []
result2 = []

for item in tqdm(data):
    input_mol = Chem.MolFromSmiles(item['input'])
    input_fp = MACCSkeys.GenMACCSKeys(input_mol)
    l1 = []
    for i in item['brics']:
        try:
            mol = Chem.MolFromSmiles(i)
            mol_fp = MACCSkeys.GenMACCSKeys(mol)
            l1.append(FingerprintSimilarity(input_fp, mol_fp))
        except:
            l1.append(0)
    l1.append(0)
    fig1.append(l1)

    l2 = []
    for i in item['ChemDual']:
        try:
            mol = Chem.MolFromSmiles(i)
            mol_fp = MACCSkeys.GenMACCSKeys(mol)
            l2.append(FingerprintSimilarity(input_fp, mol_fp))
        except:
            l2.append(0)
    fig2.append(l2)


for f1,f2 in zip(fig1, fig2):
    r1 = []
    r2 = []
    for i in step:
        r1.append(max(f1[:i]))
        r2.append(max(f2[:i]))
    result1.append(r1)
    result2.append(r2)

def calc(res):
    ans = [0, 0, 0, 0]
    for i in res:
        for j in range(len(ans)):
            ans[j] += i[j]
    for i in range(len(ans)):
        ans[i] /= len(res)
    return ans

print(calc(result1))   
print(calc(result2))            