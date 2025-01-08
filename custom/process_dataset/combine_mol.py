import os

mol_dir = "custom/mol"
output_files = "custom/mol/small_molecule.txt"
smiles_list = []
for item in os.listdir(mol_dir):
    if not item.endswith("smi"):
        continue
    for line in open(os.path.join(mol_dir, item)).readlines()[1:]:
        smiles = line.split(" ")[0]
        smiles_list.append(smiles + "\n")


with open(output_files, "w") as f:
    f.writelines(smiles_list)