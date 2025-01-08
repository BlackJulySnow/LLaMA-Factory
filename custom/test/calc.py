from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import FingerprintSimilarity

def calculate_tanimoto(smiles1, smiles2):
    # 生成分子对象
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    # 计算MACCS指纹
    fp1 = AllChem.GetMACCSKeysFingerprint(mol1)
    fp2 = AllChem.GetMACCSKeysFingerprint(mol2)
    print(fp1,fp2)
    # 计算Tanimoto相似度
    similarity = FingerprintSimilarity(fp1, fp2)
    return similarity

# 示例
smiles1 = "CCO"  # 乙醇
smiles2 = "CC(=O)O"  # 醋酸
tanimoto = calculate_tanimoto(smiles1, smiles2)
print(f"Tanimoto相似度: {tanimoto}")