import json
import pandas as pd
from sklearn.model_selection import train_test_split

def split_smiles():
    # 读取 JSON 文件
    with open('data\smiles_all.json', 'r') as file:
        data = json.load(file)

    # 将 JSON 数据转换为 DataFrame
    df = pd.DataFrame(data)
    # print(len(df))
    # # 按 8:2 的比例分割数据
    train_df, test_df = train_test_split(df, test_size=0.005, random_state=42)


    # # 将分割后的 DataFrame 转换为字典
    train_data = train_df.to_dict(orient='records')
    test_data = test_df.to_dict(orient='records')

    # # 将训练集和测试集分别写入 JSON 文件
    with open('data\smiles_train.json', 'w') as train_file:
        json.dump(train_data, train_file, indent=4)

    with open('data\smiles_test.json', 'w') as test_file:
        json.dump(test_data, test_file, indent=4)

def split_compound():
    # 读取 JSON 文件
    with open('data\compound_all.json', 'r') as file:
        data = json.load(file)

    # 将 JSON 数据转换为 DataFrame
    df = pd.DataFrame(data)
    # print(len(df))
    # # 按 8:2 的比例分割数据
    train_df, test_df = train_test_split(df, test_size=0.005, random_state=42)


    # # 将分割后的 DataFrame 转换为字典
    train_data = train_df.to_dict(orient='records')
    test_data = test_df.to_dict(orient='records')

    # # 将训练集和测试集分别写入 JSON 文件
    with open('data\compound_train.json', 'w') as train_file:
        json.dump(train_data, train_file, indent=4)

    with open('data\compound_test.json', 'w') as test_file:
        json.dump(test_data, test_file, indent=4)

split_compound()
# split_smiles()