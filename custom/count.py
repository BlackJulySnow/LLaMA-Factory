import pandas as pd
import ast
import json
import glob
dir = "custom/smiles"
instruction = "Decomposing Molecule Parts"

file_paths = glob.glob(f'{dir}/DOWN*.csv')
print(file_paths[0])
header_df = pd.read_csv(file_paths[0], low_memory=False, sep=";", on_bad_lines='warn', quotechar='"')
header = header_df.columns.tolist()

# 读取所有文件，并将header应用于其他文件
dataframes = [pd.read_csv(file, low_memory=False, sep=";", on_bad_lines='warn', quotechar='"', header=0 if i == 0 else None, names=header) for i, file in enumerate(file_paths)]
# 合并数据框
df = pd.concat(dataframes, ignore_index=True)

# 去重并保存结果
print("Total rows before drop duplicates:", len(df))
df = df.drop_duplicates(subset='Smiles')
df.to_csv(f"{dir}/smiles.csv", index=False, encoding="utf-8")
print("Total rows after drop duplicates:", len(df))

# # df = pd.read_csv(f"{dir}/smiles0.csv", low_memory=False)
# data = []
# count = 0
# for index, row in df.iterrows():
#     smiles = row['Smiles']
#     brics = row['Brics']
#     # 在这里处理每一行的数据
#     q = f'What parts can the molecule "{smiles}" be decomposed into?'
#     if pd.isna(brics):
#         continue
#     ans_list = ast.literal_eval(brics)
#     if len(ans_list) == 0:
#         continue
#     ans = f'We can decompose the molecule with the SMILES notation "{smiles}" into these parts:\n'
#     for idx,item in enumerate(ans_list):
#         ans += f"{idx + 1}. " + item + "\n"
#     data.append({
#         "instruction": instruction,
#         "input": q,
#         "output": ans.strip()
#     })
#     count += 1
# with open(f"{dir}/smiles_all.json", "w") as json_file:
#     json.dump(data, json_file, indent=4)
# print(count)