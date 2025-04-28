"""
Create a Shannon code for a given distribution of the source
alphabet.
    - generate the lengths l_i
    - create a tree with 2^l_max nodes
    - delete the ones with the lengths
"""

from math import log2, ceil
from itertools import product
from claude import shannon_fano_codebook as claude
probs = [0.25, 0.2, 0.12, 0.3,  0.08, 0.05]


def split(group):
    """
    Given a group of (symbol, probability, code) tuples, split it into 2 groups of almost
    equal probability. This is done by picking the split that results in the smallest
    difference of subgroup probabilities.
    """
    best_diff = float("inf")
    split_index = 0
    total = sum(x[1] for x in group)
    running_sum = 0

    for i in range(len(group)-1):
        running_sum += group[i][1]
        diff = abs(running_sum - (total - running_sum))

        if diff < best_diff:
            best_diff = diff
            split_index = i+1

    left = group[:split_index]
    right = group[split_index:]
    for l in left:
        l[2] += "0"
    for r in right:
        r[2] += "1"

    return left, right

def shannon_fano_codebook(probs):
    # make the probabilities into tuples for handling
    tuples = [(i, prob) for i, prob in enumerate(probs)]
    tuples = sorted(tuples, key=lambda x:x[1], reverse=True)

    groups = [[[id, prob, ""] for id, prob in tuples]]
    codebook = [None for _ in range(len(probs))]

    while len(groups) > 0:
        group = groups.pop()
        for subgroup in split(group):
            if len(subgroup) == 1:
                codebook[subgroup[0][0]] = subgroup[0][2]
            else:
                groups.insert(0, subgroup)

    return codebook

def shannon_fano_codebook_seq(n, p):
    """
    Generate a Shannon-Fano codebook for all binary sequences of length n, where 1 appears with
    probability p.
    """
    probs = []
    sequences = ["".join(x) for x in product(("0","1"), repeat=n)]
    for seq in sequences:
        ones = sum([1 for i in seq if i=="1"])
        probs.append(p**ones * (1-p)**(n-ones))
    
    return probs, shannon_fano_codebook(probs), claude(probs)



print(shannon_fano_codebook(probs))

probs, mine, claudes = shannon_fano_codebook_seq(10, 0.1)

print(len(mine), len(claudes))

errors = []
for m, c in zip(mine, claudes):
    if m != c[1]:
        errors.append((c[0], m, c[1], len(m), len(c[1])))
for error in errors:
    print(error)

print(mine[:10])
print([x[1] for x in claudes[:10]])

mine_avg = sum([len(x)*probs[i] for i, x in enumerate(mine)])
claude_avg = sum([len(x)*probs[i] for i, x in claudes])

print(f"my avg: {mine_avg}\ncl avg: {claude_avg}")
