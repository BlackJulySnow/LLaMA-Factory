from tqdm import tqdm
import json

path = "custom/dataset/ChemDual/recombination_train.json"

with open(path) as f :
    data = json.load(f)
result = []
for i in tqdm(data):
    i['instruction'] = "Help me to recombine these molecules."
    result.append(i)
with open(path, "w") as f:
    json.dump(result, f, indent=2)