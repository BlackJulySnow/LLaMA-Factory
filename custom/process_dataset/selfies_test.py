import selfies as sf
import random
from rdkit.Chem import BRICS
from rdkit import Chem

smiles = "CCC(=O)OCCOc1ccccc1"
mol = Chem.MolFromSmiles(smiles)

res = sorted(BRICS.BRICSDecompose(mol))

ms = BRICS.BRICSBuild(fragms)