from rdkit import Chem

# 含虚拟原子的 SMILES 表示
smiles_with_virtual_atoms = "N1CCc2ccc(=O)n([9*])c2CC1"

# 转换为分子对象
mol = Chem.MolFromSmiles(smiles_with_virtual_atoms)

# 重新生成合法的 SMILES
smiles_without_virtual_atoms = Chem.MolToSmiles(mol)
print(smiles_without_virtual_atoms)