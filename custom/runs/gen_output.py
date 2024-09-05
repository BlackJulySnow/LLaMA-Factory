import transformers
import torch
import json
from tqdm import tqdm
from tqdm.contrib import tzip
from datasets import load_dataset
from transformers.pipelines.pt_utils import KeyDataset

name = "mol_forward_reaction_prediction_test"
model = "llama_w_o_retrosynthesis"
model_id = f"C:\\Users\\B1GGersnow\\LLaMa-Factory\\saves\\ChemDual\\additional\\{model}\\model"
data_files = f"custom/dataset/ChemDual/{name}.json"
result_file = f"custom/result/{model}_{name}.jsonl"
batch_size = 1

pipeline = transformers.pipeline(
    "text-generation",
    model=model_id,
    model_kwargs={"torch_dtype": torch.bfloat16},
    device_map="auto",
    batch_size=batch_size,
)

# Load the JSON dataset
dataset = load_dataset("json", data_files=data_files)["train"]


def delete(i, idx):
    if idx >= 10:
        return False
    else:
        return True


# dataset = dataset.filter(delete, with_indices=True)


# Prepare prompts using the tokenizer's chat template function
def create_prompts(batch):
    prompts = []
    for instruction, input_text in zip(batch["instruction"], batch["input"]):
        messages = [
            {"role": "system", "content": instruction},
            {"role": "user", "content": input_text},
        ]
        prompt = pipeline.tokenizer.apply_chat_template(
            [messages], tokenize=False, add_generation_prompt=True
        )[0]
        prompts.append(prompt)
    return {"prompt": prompts}


dataset = dataset.map(create_prompts, batched=True, batch_size=batch_size)
terminators = [
    pipeline.tokenizer.eos_token_id,
    pipeline.tokenizer.convert_tokens_to_ids("<|eot_id|>"),
]

output = []
for out in tqdm(
    pipeline(
        KeyDataset(dataset, "prompt"),
        eos_token_id=terminators,
        max_new_tokens=512,
        do_sample=True,
        temperature=0.95,
        top_p=0.7,
        num_return_sequences=1,
    ),
    total=len(dataset),
):
    output.append(out[0]["generated_text"])
lines = []
for idx, out in enumerate(tqdm(output, desc="write results")):
    predict = out[len(dataset[idx]["prompt"]) :]
    label = dataset[idx]["output"]
    prompt = dataset[idx]["prompt"]
    lines.append(
        json.dumps({"predict": predict, "label": label, "prompt": prompt}) + "\n"
    )

with open(result_file, "w") as f:
    f.writelines(lines)
