import pandas as pd
import ast
import json
import re
from tqdm import tqdm
import sys

sys.path.append("custom")
from utils import valid_smiles

dir = "custom/smiles"


df = pd.read_csv(f"{dir}/processed_smiles.csv")
df = df.dropna(subset=["brics"])


def remove_pattern(text):
    # 使用正则表达式匹配并去除 "[数字*]" 模式
    result = re.sub(r"\[\d+\*\]", "C", text)
    return result.replace("'", '"')


def generate_brics():
    instruction = "Help me to decompose molecule into parts."
    data = []
    count = 0
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        smiles = row["Smiles"]
        brics = row["brics"]
        # 在这里处理每一行的数据
        q = smiles
        if pd.isna(brics):
            continue
        ans_list = ast.literal_eval(brics)
        if len(ans_list) == 0:
            continue
        ans = ".".join(ans_list)
        data.append({"instruction": instruction, "input": q, "output": ans.strip()})
        count += 1
    with open(f"{dir}/smiles_all.json", "w") as json_file:
        json.dump(data, json_file, indent=2)


def generate_compound():
    instruction = "Help me to synthesize molecules."
    data = []
    count = 0
    valid = 0
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        smiles = row["Smiles"]
        brics = remove_pattern(row["brics"])
        # 在这里处理每一行的数据
        if pd.isna(brics):
            continue
        q_list = ast.literal_eval(brics)
        q = ".".join(q_list)
        if valid_smiles(q):
            valid += 1
        ans = smiles
        data.append({"instruction": instruction, "input": q, "output": ans})
        count += 1
    print(valid, count)
    with open(f"{dir}/compound_all.json", "w") as json_file:
        json.dump(data, json_file, indent=2)


generate_compound()
generate_brics()
