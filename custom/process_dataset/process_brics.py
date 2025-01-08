import re
import pandas as pd
from rdkit import Chem
from custom_brics import BRICSDecompose
# from rdkit.Chem.BRICS import BRICSDecompose

# 读取CSV文件
df = pd.read_csv("custom\dataset\processed_smiles.csv")
# def remove_pattern(text):
#     # 使用正则表达式匹配 "[数字*]" 或 "([数字*])" 模式，并替换为 "C"
#     result = re.sub(r"n\[\d+\*\]", "N", text)
#     result = re.sub(r"\[\d+\*\]|\(\[\d+\*\]\)", "", text)
#     return result.replace("'", '"')
def valid_smiles(smiles):
    # 判断SMILES是否有效
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return True
        else:
            return False
    except:
        return False
# 取df前1000行数据进行测试
df = df.head(1000)

for idx, row in df.iterrows():
    smiles = row['Smiles']
    mol = Chem.MolFromSmiles(smiles)
    # fragments = BRICSDecompose(mol)
    framgents = BRICSDecompose(mol, smiles)
    
    df.at[idx, 'brics'] = '.'.join(framgents)

df.to_csv("custom\dataset\processed_brics.csv", index=False)
