import pandas as pd
import json

# 读取Parquet文件
df = pd.read_parquet("custom/train-00000-of-00001.parquet")
# 将DataFrame转换为JSON格式
json_data = df.to_json(orient='records')
data = json.loads(json_data)
# # 将JSON字符串写入文件
with open("custom/ChemData700K.json", "w") as json_file:
    json.dump(data, json_file, indent=4)

print("JSON数据已成功写入文件ChemData700K.json")
