import json

synthesize = []
decompose = []
with open(
    "saves/LLaMA3-8B/lora/eval_2024-08-09-07-19-31/generated_predictions.jsonl"
) as f:
    for line in f.readlines():
        if line.count("synthesize") == 1:
            synthesize.append(json.loads(line))
        else:
            decompose.append(json.loads(line))

with open(
    "saves/LLaMA3-8B/lora/eval_2024-08-09-07-19-31/result_synthesize.json", "w"
) as f:
    f.write(json.dumps(synthesize, indent=2))

with open(
    "saves/LLaMA3-8B/lora/eval_2024-08-09-07-19-31/result_decompose.json", "w"
) as f:
    f.write(json.dumps(decompose, indent=2))
print(len(synthesize), len(decompose))