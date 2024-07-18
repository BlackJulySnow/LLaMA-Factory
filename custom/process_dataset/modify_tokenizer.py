import json

with open("C:\\Users\\B1GGersnow\\Desktop\\model\\LLama3-8B-Instruct\\tokenizer.json.bak",'r', encoding='utf-8') as f1:
    with open("C:\\Users\\B1GGersnow\\Desktop\\model\\LLama3-8B-Instruct\\vocab_bpe_300_result.txt") as f2:
        res = json.loads(f1.read())
        for idx, item in enumerate(f2.readlines()):
            res['added_tokens'][idx]['content'] = item[:-1]
        with open("C:\\Users\\B1GGersnow\\Desktop\\model\\LLama3-8B-Instruct\\tokenizer_smiles.json", 'w', encoding='utf-8') as f3:
            json.dump(res, f3, indent=2, ensure_ascii=False)
