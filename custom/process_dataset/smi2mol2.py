import json
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import BRICS, AllChem
import random
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import FingerprintSimilarity

with open("custom\dataset\ChemDual_fragment_test.json") as f:
    data = json.load(f)

for idx,item in tqdm(enumerate(data[:100])):
    l = item['output'].split('.')
    random.seed(127)
    fragms = [Chem.MolFromSmiles(x) for x in l]
    random.seed(0xf00d)
    ms = BRICS.BRICSBuild(fragms)
    input_mol = Chem.MolFromSmiles(item['input'])
    Chem.SanitizeMol(input_mol)
    item['input'] = Chem.MolToSmiles(input_mol)
    # input_fp = MACCSkeys.GenMACCSKeys(input_mol)
    # max_similarity = -1
    # best_mol = None
    l = []
    try:
        for i in range(10):
            # 获取下一个重组分子
            mol = next(ms)
            if mol is not None:
                Chem.SanitizeMol(mol)
                # mol_fp = MACCSkeys.GenMACCSKeys(mol)
                # 计算与原始分子的相似度
                # similarity = FingerprintSimilarity(input_fp, mol_fp)
                l.append(Chem.MolToSmiles(mol))
                # 更新最大相似度和最佳分子
                # if similarity > max_similarity:
                    # max_similarity = similarity
                    # best_mol = mol
            else:
                break
    except StopIteration:
        pass
        # 迭代器耗尽时会引发 StopIteration
    # if best_mol:
    #     best_smiles = Chem.MolToSmiles(best_mol)
    #     item['brics_recombine'] = best_smiles
    #     item['similarity'] = max_similarity
    #     # print(f"Best reassembled molecule for {idx + 1} with max similarity {max_similarity}: {best_smiles}")
    # else:
    #     print(f"No reassembled molecules found for {idx + 1}")
    item['brics'] = l

with open("custom\dataset\ChemDual_fragment_test_brics.json", "w") as f:
    json.dump(data, f, indent=2)