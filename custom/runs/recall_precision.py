from collections import Counter
from rdkit import Chem


def tokenize(text):
    return text.split()


def lcs_length(x, y):
    m = len(x)
    n = len(y)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m):
        for j in range(n):
            if x[i] == y[j]:
                dp[i + 1][j + 1] = dp[i][j] + 1
            else:
                dp[i + 1][j + 1] = max(dp[i + 1][j], dp[i][j + 1])

    return dp[m][n]


def calculate_memgp(candidate, reference, b = 1):
    candidate_tokens = tokenize(candidate)
    reference_tokens = tokenize(reference)
    rd_error = [not rd_exist(i) for i in reference_tokens]
    candidate_counts = Counter(candidate_tokens)
    reference_counts = Counter(reference_tokens)
    matches = sum((candidate_counts & reference_counts).values())
    bp = (1 - sum(rd_error) / len(reference_tokens)) if len(reference_tokens) > 0 else 0
    lcs_score = 0
    for i in reference_tokens:
        res = 0
        for j in candidate_tokens:
            res = max(res, lcs_length(i,j))
        lcs_score += res
    total_gen_len = sum(len(i) for i in reference_tokens)
    lcs_score = lcs_score / total_gen_len if total_gen_len > 0 else 0
    precision = matches / len(candidate_tokens) if len(candidate_tokens) > 0 else 0
    recall = matches / len(reference_tokens) if len(reference_tokens) > 0 else 0
    if precision + recall > 0:
        fmeasure = (b * b + 1) * precision * recall / (b * b * precision + recall)
    else:
        fmeasure = 0
    return {"precision": precision, "recall": recall, "fmeasure": fmeasure, "memgp": bp * (fmeasure * 0.5 + lcs_score * 0.5)}


def calculate_rouge_l(candidate, reference):
    candidate_tokens = tokenize(candidate)
    reference_tokens = tokenize(reference)

    lcs = lcs_length(candidate_tokens, reference_tokens)

    precision = lcs / len(candidate_tokens) if len(candidate_tokens) > 0 else 0
    recall = lcs / len(reference_tokens) if len(reference_tokens) > 0 else 0
    if precision + recall > 0:
        f1_score = 2 * precision * recall / (precision + recall)
    else:
        f1_score = 0

    return {"precision": precision, "recall": recall, "f1_score": f1_score}


def rd_exist(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return not mol is None


if __name__ == "__main__":
    # 示例句子
    candidate = "[3*]O[3*] [5*]N[5*] [1*]C([6*])=O [4*]C[8*] [16*]c1ccccc1 [16*]c1cccc(Cl)c1C [14*]c1nc([14*])c(C)o1 [16*]c1ccc([16*])cc1"
    reference = "[3*]O[3*] [5*]N[5*] [1*]C([6*])=O [4*]C[8*] [16*]c1ccccc1 [16*]c1cccc(Cl)c1C [14*]c1nc([14*])c(C)o1 c"

    # 计算1-gram的精度和召回率
    precision, recall, fmeasure, memgp = calculate_memgp(
        candidate, reference
    ).values()
    print(f"recision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F-measure Score: {fmeasure:.4f}")
    print(f"memgp Score: {memgp:.4f}")

    # 计算ROUGE-L
    precision_rouge_l, recall_rouge_l, f1_rouge_l = calculate_rouge_l(
        candidate, reference
    ).values()
    print(f"ROUGE-L Precision: {precision_rouge_l:.4f}")
    print(f"ROUGE-L Recall: {recall_rouge_l:.4f}")
    print(f"ROUGE-L F1 Score: {f1_rouge_l:.4f}")
