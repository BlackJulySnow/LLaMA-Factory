import json
from tqdm import tqdm
small_molecules = []
with open("custom/mol/small_molecule.txt", "r") as f:
    for line in f.readlines():
        small_molecules.append(line[:-1])
tot = 0
with open("custom/dataset/ChemDual_fragment_train.json" ,'r') as f:
    data = json.load(f)
    for item in tqdm(data):
        if item['output'] in small_molecules:
            tot += 1
            print(item['output'])
print(f"{tot}/{len(data)}")
