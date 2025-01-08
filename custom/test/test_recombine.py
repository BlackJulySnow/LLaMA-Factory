from rdkit import Chem
from rdkit.Chem import AllChem

# 示例：SMILES 列表
smiles_list = [
    "CCO",  # 乙醇
    "CCOCC",  # 乙基乙醇
    "CCC"  # 丙烷
]

# 创建一个空的分子列表
mols = []

# 遍历 SMILES 列表，生成分子对象
for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)  # 从 SMILES 创建分子
    if mol is not None:
        AllChem.Compute2DCoords(mol)  # 生成 2D 坐标
        mols.append(mol)

# 将分子保存为 SDF 格式
w = Chem.SDWriter('output.sdf')  # 创建 SDF 写入对象
for mol in mols:
    w.write(mol)  # 写入每个分子到 SDF 文件中
w.close()  # 关闭写入对象
