import json
import re
from tqdm import tqdm

with open("custom/dataset/smiles_test.json") as f:
    dataset = json.load(f)
    dataset = dataset[:1000]

with open("custom/dataset/result_brics.json", "r") as f:
    results = json.load(f)


def remove_pre(s):
    pattern = re.compile(r"^\d+\.")
    return pattern.sub("", s).lstrip().strip()


def remove_pattern(text):
    # 使用正则表达式匹配并去除 "[数字*]" 模式
    result = re.sub(r"\[\d+\*\]", "", text)
    return result.replace("'", '"')


origin = []
for item in dataset:
    l = item["output"].split("\n")[1:]
    res = [remove_pattern(remove_pre(i)) for i in l]
    res = sorted(res)
    origin.append(res)

output = []
for item in results:
    l = item[0].split("\n")
    pattern = re.compile(r"^\d+\..*")
    filtered_list = [remove_pattern(remove_pre(s)) for s in l if pattern.match(s)]
    filtered_list = sorted(filtered_list)
    output.append(filtered_list)

from rouge_score import rouge_scorer

scorer = rouge_scorer.RougeScorer(["rouge1", "rougeL"], use_stemmer=True)
total_scores = {"rouge1": 0, "rougeL": 0}
num_pairs = len(origin)

for i, j in tqdm(zip(origin, output)):
    str_origin = " ".join(i)
    str_output = " ".join(j)
    scores = scorer.score(str_origin, str_output)
    for key in total_scores:
        total_scores[key] += scores[key].fmeasure

average_scores = {key: value / num_pairs for key, value in total_scores.items()}

print(average_scores)
