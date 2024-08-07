import json
from valid_smiles import valid_smiles
import csv
from nltk.translate.bleu_score import sentence_bleu, SmoothingFunction

with open("custom/dataset/compound_test.json") as f:
    dataset = json.load(f)
    dataset = dataset[:1000]

with open("custom/dataset/result_compound.json", "r") as f:
    output = json.load(f)

origin = []
for item in dataset:
    head = "If you have these small molecules, you can synthesize "
    l = item["output"][len(head) :]
    origin.append(l)

# with open("caption2smiles_example.txt", "w", newline="") as csvfile:
#     csvwriter = csv.writer(csvfile, delimiter="\t")
#     # 写入表头
#     csvwriter.writerow(["description", "ground_smiles", "output_smiles"])
#     # 写入每一行的数据
#     for i, j in zip(origin, output):
#         csvwriter.writerow(["description", i, j[0]])

error_item = []
def calc(accuracy_k: int):
    accuracy_count = 0
    valid_count = 0
    for item_origin, item_output in zip(origin, output):
        if sum(valid_smiles(item_output[:accuracy_k])) == 0 and item_origin not in error_item:
            error_item.append(item_origin)
        valid_count += sum(valid_smiles(item_output[:accuracy_k]))
        accuracy_count += item_origin in item_output[:accuracy_k]
    return accuracy_count, valid_count


accuracy_count1, valid_count1= calc(1)
accuracy_count3, valid_count3= calc(3)
accuracy_count5, valid_count5= calc(5)
print(error_item)
print(f"accuracy1: {accuracy_count1 / 10:.2f}%, validity1: {valid_count1 / 10 / 1:.2f}%")
print(f"accuracy3: {accuracy_count3 / 10:.2f}%, validity3: {valid_count3 / 10 / 3:.2f}%")
print(f"accuracy5: {accuracy_count5 / 10:.2f}%, validity5: {valid_count5 / 10 / 5:.2f}%")
