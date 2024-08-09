from rdkit import Chem
from rdkit.Chem import Draw


def draw(smiles):
    # 解析SMILES字符串为分子对象
    molecule = Chem.MolFromSmiles(smiles)

    # 检查解析是否成功
    if molecule is None:
        print("无法解析SMILES字符串。")
    else:
        # 可视化分子
        img = Draw.MolToImage(molecule)
        img.show()


def valid_smiles(smiles: str) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return True
        else:
            return False
    except:
        return False


if __name__ == "__main__":
    # Example usage
    smiles_list = "C1CCCCC1"
    print(valid_smiles(smiles_list))
