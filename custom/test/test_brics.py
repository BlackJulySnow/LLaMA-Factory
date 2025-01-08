import re
from rdkit.Chem import BRICS
from rdkit import Chem
def remove_pattern(text):
    # 使用正则表达式匹配并去除 "[数字*]" 模式
    result = re.sub(r"\[\d+\*\]", "C", text)
    return result.replace("'", '"')
mol = Chem.MolFromSmiles("NC(=O)CN1CCC(N)CC1")

res = sorted(BRICS.BRICSDecompose(mol))
print('.'.join(res))
print(remove_pattern('.'.join(res)))
# res = ['CC(C)=O', 'C[C@@H]1C[C@H]2C[C@@H]([15*])[C@H]2C1', '[4*]CCCC', '[5*]N[5*]', '[6*]C(=O)O']
# fragms = [Chem.MolFromSmiles(x) for x in sorted(res)]
# # print(len(res))
# # print(".".join(res))
# ms = BRICS.BRICSBuild(fragms)
# prod = [next(ms) for x in range(3)]
# for i in prod:
#     print(Chem.MolToSmiles(i))

