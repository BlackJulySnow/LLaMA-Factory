from rdkit import Chem
from rdkit.Chem import Draw

# SMILES字符串
smiles = '[C][C][O][C][=Branch1][C][=O][C][NH1][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=Ring1][=Branch2][C][C][C][Ring1][=Branch1][O]'

# 解析SMILES字符串为分子对象
molecule = Chem.MolFromSmiles(smiles)

# 检查解析是否成功
if molecule is None:
    print("无法解析SMILES字符串。")
else:
    # 可视化分子
    img = Draw.MolToImage(molecule)
    img.show()
