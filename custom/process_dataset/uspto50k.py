import pandas as pd

# 设置Parquet文件路径
parquet_file_path = 'custom/dataset/uspto-50k.parquet'

# 读取Parquet文件
df = pd.read_parquet(parquet_file_path)

# 打印数据集信息
print("DataFrame:")
print(df.head())

# 打印数据集统计信息
print("\nDataFrame Info:")
print(df.info())

# 打印数据集描述统计
print("\nDataFrame Description:")
print(df.describe())
