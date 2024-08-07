import argparse
import csv

import os.path as osp

import numpy as np

from transformers import AutoTokenizer

from nltk.translate.bleu_score import corpus_bleu
from nltk.translate.meteor_score import meteor_score
from rouge_score import rouge_scorer

def evaluate(text_model, input_file, text_trunc_length):
    outputs = []

    with open(osp.join(input_file), encoding='utf8') as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for n, line in enumerate(reader):
            outputs.append((line['ground_smiles'], line['output_smiles']))

    text_tokenizer = AutoTokenizer.from_pretrained(text_model)

    bleu_scores = []
    meteor_scores = []

    references = []
    hypotheses = []

    for i, (gt, out) in enumerate(outputs):

        if i % 100 == 0: print(i, 'processed.')


        gt_tokens = text_tokenizer.tokenize(gt, truncation=True, max_length=text_trunc_length,
                                            padding='max_length')
        gt_tokens = list(filter(('<|eot_id|>').__ne__, gt_tokens))

        out_tokens = text_tokenizer.tokenize(out, truncation=True, max_length=text_trunc_length,
                                            padding='max_length')
        out_tokens = list(filter(('<|eot_id|>').__ne__, out_tokens))

        references.append([gt_tokens])
        hypotheses.append(out_tokens)

        mscore = meteor_score([gt_tokens], out_tokens)
        meteor_scores.append(mscore)

    bleu2 = corpus_bleu(references, hypotheses, weights=(.5,.5))
    bleu4 = corpus_bleu(references, hypotheses, weights=(.25,.25,.25,.25))

    print('BLEU-2 score:', bleu2)
    print('BLEU-4 score:', bleu4)
    _meteor_score = np.mean(meteor_scores)
    print('Average Meteor score:', _meteor_score)

    scorer = rouge_scorer.RougeScorer(['rouge1', 'rouge2', 'rougeL'])

    rouge_scores = []

    references = []
    hypotheses = []

    for i, (gt, out) in enumerate(outputs):

        rs = scorer.score(out, gt)
        rouge_scores.append(rs)

    print('ROUGE score:')
    rouge_1 = np.mean([rs['rouge1'].fmeasure for rs in rouge_scores])
    rouge_2 = np.mean([rs['rouge2'].fmeasure for rs in rouge_scores])
    rouge_l = np.mean([rs['rougeL'].fmeasure for rs in rouge_scores])
    print('rouge1:', rouge_1)
    print('rouge2:', rouge_2)
    print('rougeL:', rouge_l)
    return bleu2, bleu4, rouge_1, rouge_2, rouge_l, _meteor_score

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--text_model', type=str, default='C:\\Users\\B1GGersnow\\Desktop\\model\\smilesv3', help='Desired language model tokenizer.')
    parser.add_argument('--input_file', type=str, default='caption2smiles_example.txt', help='path where test generations are saved')
    parser.add_argument('--text_trunc_length', type=str, default=512, help='tokenizer maximum length')
    args = parser.parse_args()
    evaluate(args.text_model, args.input_file, args.text_trunc_length)