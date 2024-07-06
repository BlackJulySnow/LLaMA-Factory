from nltk.translate.bleu_score import sentence_bleu
from rouge_score import rouge_scorer

def calculate_rouge_scores(generated, reference):
    scorer = rouge_scorer.RougeScorer(['rouge1', 'rouge2', 'rougeL'], use_stemmer=True)
    scores = scorer.score(' '.join(reference), ' '.join(generated))
    return scores

# 示例序列
generated_sequence = ['CC(C)C', 'CC(C)C', 'CC(C)C']
reference_sequence = ['CC(C)C', 'CC(C)C', 'CC(C)C', 'c1ccccc1']

# 评价
results = calculate_rouge_scores(generated_sequence, reference_sequence)
print(results)
