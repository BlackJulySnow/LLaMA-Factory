import json

path = "custom/dataset/ChemDual/mol_forward_reaction_prediction_train.json"

with open(path) as f :
    data = json.load(f)
instructions = set()
for i in data:
    instructions.add(i['instruction'])
print(list(instructions))