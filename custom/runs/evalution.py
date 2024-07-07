import transformers
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM

# 加载分词器


# 加载模型
# model_id = "C:\\Users\\B1GGersnow\\Desktop\\test"
model_name = "C:\\Users\\B1GGersnow\\Desktop\\smilesv1"

# print(tokenizer.eos_token_id)
pipeline = transformers.pipeline(
    "text-generation",
    model=model_name,
    model_kwargs={"torch_dtype": torch.bfloat16},
    device="cuda",
)
 
instruction = "Decomposing Molecule Parts"
user_input = "What parts can the molecule \"CC[C@H]1OC(=O)[C@H](C)C(=O)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(OC)C[C@@H](C)/C(=N/OCOCCOC)[C@H](C)[C@@H](O)[C@]1(C)O\" be decomposed into?"
 
messages = [
    {"role": "system", "content": instruction},
    {"role": "user", "content": user_input},
]
 
prompt = pipeline.tokenizer.apply_chat_template(
    messages,
    tokenize=False,
    add_generation_prompt=True
)
 
 
outputs = pipeline(
    prompt,
    max_new_tokens=256,
    do_sample=True,
    temperature=0.6,
    top_p=0.9,
)
print(outputs[0]["generated_text"][len(prompt):])
