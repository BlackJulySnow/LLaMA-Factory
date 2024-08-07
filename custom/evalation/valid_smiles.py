from rdkit import Chem

def valid_smiles(smiles_list: list[str]) -> list[bool]:
    valid_list = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                valid_list.append(True)
            else:
                valid_list.append(False)
        except:
            valid_list.append(False)
    return valid_list

if __name__ == "__main__":
    # Example usage
    smiles_list = ["CCO", "invalid_smiles", "C1CCCCC1", "O=C=O"]
    print(valid_smiles(smiles_list))
