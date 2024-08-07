import transformers
import torch
import json
from tqdm import tqdm

model_id = "C:\\Users\\B1GGersnow\\Desktop\\model\\smilesv3"
pipeline = transformers.pipeline(
    "text-generation",
    model=model_id,
    model_kwargs={"torch_dtype": torch.bfloat16},
    device_map="auto",
)

# 加载自定义数据集
with open("custom/dataset/compound_test.json", "r") as file:
    dataset = json.load(file)
dataset = dataset[:1000]


def calc_sim(output, original):
    from rdkit import Chem
    from collections import Counter

    for i in output:
        try:
            mol = Chem.MolFromSmiles(i)
            atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
            atom_count = Counter(atom_symbols)
            print(atom_count)
        except:
            print("error")
            pass


def generate_text(
    instruction,
    user_input,
    max_new_tokens=256,
    do_sample=True,
    temperature=1.5,
    top_p=0.4,
    num_return_sequences=3,
):
    messages = [
        {"role": "system", "content": instruction},
        {"role": "user", "content": user_input},
    ]
    prompt = pipeline.tokenizer.apply_chat_template(
        messages, tokenize=False, add_generation_prompt=True
    )
    outputs = pipeline(
        prompt,
        max_new_tokens=max_new_tokens,
        do_sample=do_sample,
        temperature=temperature,
        top_p=top_p,
        num_return_sequences=num_return_sequences,
    )
    return [output["generated_text"][len(prompt)] for output in outputs]


num_return_sequences = 5
res1 = []

for item in tqdm(dataset, desc="Processing items"):
    instruction = item["instruction"]
    user_input = item["input"]
    expected_output = item["output"]
    generated_texts = generate_text(
        instruction,
        user_input,
        num_return_sequences=num_return_sequences,
        temperature=0.6,
        top_p=0.9,
    )
    res1.append(generated_texts)
with open("custom/dataset/result_compound.json", "w") as result_file:
    json.dump(res1, result_file)
    # for i, generated_text in enumerate(generated_texts, 1):
    #     print(f"Generated {i}: {generated_text}")
    # print(f"Expected: {expected_output}\n")
