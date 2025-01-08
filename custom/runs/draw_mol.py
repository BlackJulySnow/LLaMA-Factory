from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

# 分子SMILES列表
smiles_list = [
    "COc1cc(NC(=O)[C@@H](O)[C@H](N)CC2CCCCC2)ccc1NC(=O)/C=C/c1cc(OC)c(OC)c(OC)c1",
    "Cc1cc(/C=C/C(=O)c2ccc(C)c(C)c2)cc(C)c1C[C@@H](N)[C@H](O)C(=O)C1CCCCC1",
    "CC[C@@H](N)[C@H](O)C(=O)c1cc(/C=C/C(=O)c2ccc(C)c(C)c2)cc(C2CCCCC2)c1C",
    "Cc1ccc(C[C@@H](N)[C@H](O)C(=O)c2cc(C)c(CC3CCCCC3)c(C)c2)cc1C(=O)CC1CCCCC1",
    "CC[C@@H](N)[C@H](O)C(=O)c1ccc(/C=C/C(=O)c2cc(C)c(C)c(C)c2)cc1C1CCCCC1",
    "CC[C@@H](N)[C@H](O)C(=O)c1ccc(C(=O)/C=C/c2cc(C)c(C)c(C)c2)cc1C1CCCCC1",
    "CC[C@@H](N)[C@H](O)C(=O)c1cc(/C=C/C(=O)c2ccc(C)c(C)c2)cc(C)c1C1CCCCC1"
]

# 创建分子对象
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# 创建图形
fig = plt.figure(figsize=(15, 12))
gs = plt.GridSpec(3, 3, figure=fig)

# 第一行的大图
ax1 = fig.add_subplot(gs[0, :])
img1 = Draw.MolToImage(mols[0])
ax1.imshow(img1)
ax1.axis('off')
ax1.set_title('a', y=-0.15)

# 第二行的三个图
for i in range(3):
    ax = fig.add_subplot(gs[1, i])
    img = Draw.MolToImage(mols[i+1], (800, 800))
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(f'{chr(97+i)}', y=-0.15)

# 第三行的三个图
for i in range(3):
    ax = fig.add_subplot(gs[2, i])
    img = Draw.MolToImage(mols[i+4], (800, 800))
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(f'{chr(97+i)}', y=-0.15)

# 调整子图间距
plt.tight_layout()

# 显示图形
plt.show()