import multiprocessing
import pandas as pd
from custom_brics import BRICSDecompose
from rdkit import Chem
df = pd.read_csv("custom/smiles/smiles.csv", low_memory=False)
# 仅取前1000条数据进行测试
df = df.head(1000).copy()
df['index'] = df.index
process_size = 16
df["bin"] = pd.qcut(df["index"], q=process_size, labels=False)


# 处理Smiles的函数
def process_smiles(bin_df, process, return_dict):
    print(f"start{process}:{len(bin_df)}")
    results = []
    for row in bin_df.itertuples(index=False):
        try:
            mol = Chem.MolFromSmiles(row.Smiles)
            res = BRICSDecompose(mol, row.Smiles)
            results.append(res)
        except Exception as e:
            print(f"error: {e}")
            
    print(f"success:{process}")
    return_dict[process] = results

if __name__ == "__main__":
    # 使用多进程处理每个Smiles
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    processes = []
    for index in range(process_size):
        bin_df = df[df["bin"] == index].copy()
        p = multiprocessing.Process(target=process_smiles, args=(bin_df, index, return_dict))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    df["brics"] = [None] * len(df)
    for key, results in return_dict.items():
        bin_indices = df[df["bin"] == key].index
        for i, idx in enumerate(bin_indices):
            df.at[idx, "brics"] = results[i]

    df.to_csv("custom/smiles/processed_smiles.csv", index=False, encoding="utf-8")
    print("Processing and insertion completed.")