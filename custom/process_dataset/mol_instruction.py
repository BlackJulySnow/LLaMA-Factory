import json
import selfies as sf
from tqdm import tqdm


def sf_encode(selfies):
    try:
        smiles = sf.decoder(selfies)
        return smiles
    except Exception:
        return None


base_path = "custom/dataset/mol_instruction/"
file_path = "reagent_prediction.json"
# 读取JSON文件
with open(base_path + file_path, "r") as file:
    data = json.load(file)

# 筛选出所有split为'test'的数据
train_data = []
test_data = []
for item in tqdm(data):
    item["input"] = sf_encode(item["input"])
    item["output"] = sf_encode(item["output"])
    if item["metadata"]["split"] == "test":
        item.pop('metadata')
        test_data.append(item)
        continue
    if item["metadata"]["split"] == "train":
        item.pop('metadata')
        train_data.append(item)
        continue


# with open(f"custom/dataset/mol_{file_path[:-5]}_train.json", "w") as f:
#     f.write(json.dumps(train_data, indent=2))
# with open(f"custom/dataset/mol_{file_path[:-5]}_test.json", "w") as f:
#     f.write(json.dumps(test_data, indent=2))
print(len(train_data), len(test_data))