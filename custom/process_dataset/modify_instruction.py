import re
from tqdm import tqdm
import json
import random

path = "custom\dataset\ChemDual_fragment_test_brics.json"

# instructions = [
#     "Please assist me in breaking down the molecule into its components.",
#     "Help me to decompose molecule into parts.",
#     "Help me to break the molecule into its individual parts.",
#     "Could you help me separate the molecule into its constituent parts?",
#     "Please help me disassemble the molecule into smaller units."
# ]
# instructions = [
#     "Help me to recombine these molecules.",
#     "Please assist me in recombining these molecules.",
#     "Help me to combine these molecules together again.",
#     "I need help recombining these molecules into a single structure.",
#     "Please help me reassemble these molecules."
# ]
def remove_pattern(text):
    # 使用正则表达式匹配并去除 "[数字*]" 模式
    result = re.sub(r"\[\d+\*\]", "C", text)
    return result.replace("'", '"')
with open(path) as f :
    data = json.load(f)
result = []
for i in tqdm(data):
    i['component'] = remove_pattern(i['output'])
    # i['instruction'] = random.choice(instructions)
    result.append(i)
with open(path, "w") as f:
    json.dump(result, f, indent=2)