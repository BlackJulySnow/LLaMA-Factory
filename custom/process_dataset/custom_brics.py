from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
import re

environs = {
  'L1': '[C;D3]([#0,#6,#7,#8])(=O)',
  'L3': '[O;D2]-;!@[#0,#6,#1]',
  'L4': '[C;!D1;!$(C=*)]-;!@[#6]',
  'L5': '[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]',
  'L6': '[C;D3;!R](=O)-;!@[#0,#6,#7,#8]',
  'L7a': '[C;D2,D3]-[#6]',
  'L7b': '[C;D2,D3]-[#6]',
  '#L8': '[C;!R;!D1]-;!@[#6]',
  'L8': '[C;!R;!D1;!$(C!-*)]',
  'L9': '[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]',
  'L10': '[N;R;$(N(@C(=O))@[C,N,O,S])]',
  'L11': '[S;D2](-;!@[#0,#6])',
  'L12': '[S;D4]([#6,#0])(=O)(=O)',
  'L13': '[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]',
  'L14': '[c;$(c(:[c,n,o,s]):[n,o,s])]',
  'L14b': '[c;$(c(:[c,n,o,s]):[n,o,s])]',
  'L15': '[C;$(C(-;@C)-;@C)]',
  'L16': '[c;$(c(:c):c)]',
  'L16b': '[c;$(c(:c):c)]',
}
reactionDefs = (
  # L1
  [
    ('1', '3', '-'),
    ('1', '5', '-'),
    ('1', '10', '-'),
  ],

  # L3
  [
    ('3', '4', '-'),
    ('3', '13', '-'),
    ('3', '14', '-'),
    ('3', '15', '-'),
    ('3', '16', '-'),
  ],

  # L4
  [
    ('4', '5', '-'),
    ('4', '11', '-'),
  ],

  # L5
  [
    ('5', '12', '-'),
    ('5', '14', '-'),
    ('5', '16', '-'),
    ('5', '13', '-'),
    ('5', '15', '-'),
  ],

  # L6
  [
    ('6', '13', '-'),
    ('6', '14', '-'),
    ('6', '15', '-'),
    ('6', '16', '-'),
  ],

  # L7
  [
    ('7a', '7b', '='),
  ],

  # L8
  [
    ('8', '9', '-'),
    ('8', '10', '-'),
    ('8', '13', '-'),
    ('8', '14', '-'),
    ('8', '15', '-'),
    ('8', '16', '-'),
  ],

  # L9
  [
    ('9', '13', '-'),  # not in original paper
    ('9', '14', '-'),  # not in original paper
    ('9', '15', '-'),
    ('9', '16', '-'),
  ],

  # L10
  [
    ('10', '13', '-'),
    ('10', '14', '-'),
    ('10', '15', '-'),
    ('10', '16', '-'),
  ],

  # L11
  [
    ('11', '13', '-'),
    ('11', '14', '-'),
    ('11', '15', '-'),
    ('11', '16', '-'),
  ],

  # L12
  # none left

  # L13
  [
    ('13', '14', '-'),
    ('13', '15', '-'),
    ('13', '16', '-'),
  ],

  # L14
  [
    ('14', '14', '-'),  # not in original paper
    ('14', '15', '-'),
    ('14', '16', '-'),
  ],

  # L15
  [
    ('15', '16', '-'),
  ],

  # L16
  [
    ('16', '16', '-'),  # not in original paper
  ],
)

smartsGps = reactionDefs
for gp in smartsGps:
    for j, defn in enumerate(gp):
        g1, g2, bnd = defn
        r1 = environs['L' + g1]
        r2 = environs['L' + g2]
        g1 = re.sub('[a-z,A-Z]', '', g1)
        g2 = re.sub('[a-z,A-Z]', '', g2)
        sma = '[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]' % (r1, bnd, r2, g1, g2)
        gp[j] = sma
reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])

def BRICSDecompose(mol, mSmi, allNodes=None, minFragmentSize=1, onlyUseReactions=None, silent=True,
                   keepNonLeafNodes=False, singlePass=False, returnMols=False, maxIterations=1000):

  if allNodes is None:
    allNodes = set()

  if mSmi in allNodes:
    return set()

  activePool = {mSmi: mol}
  allNodes.add(mSmi)
  foundMols = {mSmi: mol}
  iteration_count = 0  # Add an iteration counter

  for gpIdx, reactionGp in enumerate(reactions):
    newPool = {}
    while activePool:
      if iteration_count >= maxIterations:  # Check if the maximum number of iterations is reached
        return None
      iteration_count += 1  # Increment the iteration counter
      matched = False
      nSmi = next(iter(activePool))
      mol = activePool.pop(nSmi)
      for rxnIdx, reaction in enumerate(reactionGp):
        if onlyUseReactions and (gpIdx, rxnIdx) not in onlyUseReactions:
          continue
        if not silent:
          print('--------')
          print(smartsGps[gpIdx][rxnIdx])
        ps = reaction.RunReactants((mol, ))
        if ps:
          if not silent:
            print(nSmi, '->', len(ps), 'products')
          for prodSeq in ps:
            seqOk = True
            # we want to disqualify small fragments, so sort the product sequence by size
            tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx) for idx, prod in enumerate(prodSeq)]
            tSeq.sort()
            for nats, idx in tSeq:
              prod = prodSeq[idx]
              try:
                Chem.SanitizeMol(prod)
              except Exception:
                continue
              pSmi = Chem.MolToSmiles(prod, 1)
              if minFragmentSize > 0:
                nDummies = pSmi.count('*')
                if nats - nDummies < minFragmentSize:
                  seqOk = False
                  break
              prod.pSmi = pSmi
            ts = [(x, prodSeq[y]) for x, y in tSeq]
            prodSeq = ts
            if seqOk:
              matched = True
              for nats, prod in prodSeq:
                pSmi = prod.pSmi
                # print('\t',nats,pSmi)
                if pSmi not in allNodes:
                  if not singlePass:
                    activePool[pSmi] = prod
                  allNodes.add(pSmi)
                  foundMols[pSmi] = prod
      if singlePass or keepNonLeafNodes or not matched:
        newPool[nSmi] = mol
    activePool = newPool
  if not (singlePass or keepNonLeafNodes):
    if not returnMols:
      res = set(activePool.keys())
    else:
      res = activePool.values()
  else:
    if not returnMols:
      res = allNodes
    else:
      res = foundMols.values()
  return res

if __name__ == "__main__":
    from rdkit import Chem
    smiles = 'CC[C@]1(O)CC[C@H]2[C@@H]3CCC4=CCCC[C@@H]4[C@H]3CC[C@@]21C'
    mol = Chem.MolFromSmiles(smiles)
    fragments = BRICSDecompose(mol, smiles)

    sorted(fragments)
    # 输出分解后的片段
    for frag in list(fragments):
        print(frag)