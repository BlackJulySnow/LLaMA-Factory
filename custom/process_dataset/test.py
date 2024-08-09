import sys
from rdkit import Chem
sys.path.append('custom')
from evalation.valid_smiles import valid_smiles
import selfies as sf
def sf_encode(selfies):
    try:
        smiles = sf.decoder(selfies)
        return smiles
    except Exception:
        return None


print(valid_smiles("[16*]c1ccc(C)cc1.[15*][C@@H]1CC=CC(=O)[C@@]2(C1=O)C(=O)N2N=C(C)C"))
# print(sf.encoder("[1*]C([6*])=O"))
# print(sf_encode("[C][C][Branch1][C][C][=O]"))
# print(sf_encode("[1*]C([6*])=O"))