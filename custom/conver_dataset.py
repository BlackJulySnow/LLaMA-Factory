import pandas as pd
import ast
import json
import re
dir = "custom/smiles"


df = pd.read_csv(f"{dir}/processed_smiles.csv")
df = df.dropna(subset=['brics'])


def remove_pattern(text):
    # 使用正则表达式匹配并去除 "[数字*]" 模式
    result = re.sub(r'\[\d+\*\]', '', text)
    return result.replace("'", '"')

def generate_brics():
    instruction = "Decomposing Molecule Parts"
    data = []
    count = 0
    for index, row in df.iterrows():
        smiles = row['Smiles']
        brics = row['brics']
        # 在这里处理每一行的数据
        q = f'What parts can the molecule "{smiles}" be decomposed into?'
        if pd.isna(brics):
            continue
        ans_list = ast.literal_eval(brics)
        if len(ans_list) == 0:
            continue
        ans = f'We can decompose the molecule with the SMILES notation "{smiles}" into these parts:'
        for idx,item in enumerate(ans_list):
            ans += f"\n{idx + 1}. {item}"
        data.append({
            "instruction": instruction,
            "input": q,
            "output": ans.strip()
        })
        count += 1
    with open(f"{dir}/smiles_all.json", "w") as json_file:
        json.dump(data, json_file, indent=4)
    print(count)

def generate_compound():
    instruction = "Help me synthesize molecules."
    data = []
    count = 0
    for index, row in df.iterrows():
        smiles = row['Smiles']
        brics = remove_pattern(row['brics'])[1:-1]
        # 在这里处理每一行的数据
        q = f'If I have these small molecules [{brics}], what compounds can I combine?'
        ans = f'If you have these small molecules, you can synthesize {smiles}'
        data.append({
            "instruction": instruction,
            "input": q,
            "output": ans.strip()
        })
        count += 1
    with open(f"{dir}/compound_all.json", "w") as json_file:
        json.dump(data, json_file, indent=4)
    print(count)

generate_compound()