"""
Typical set stuff.
"""
import random 
from math import log2
from itertools import product

# good example with p=0.3, n=15 and EPSILON=0.25 

EPSILON = 0.25

def ones(s: str):
    return sum([1 for i in s if i=="1"])

def probability(sequence: str):
    _ones = ones(sequence)
    zeros = len(sequence)-_ones
    return p**_ones * (1-p)**zeros

def binary_source_entropy(p):
    return -p*log2(p) - (1-p)*log2(1-p) 

p = 0.3     # the more even p is, the higher the acceptable epsilon for a fixed n (and the working n for a fixed epsilon)
n = 15      # 20->2sec, 23->20sec, 24->44sec, 25 apagoreutiko

print(f"assuming a source with p(1) = {p} (entropy of {binary_source_entropy(p)}), sequences of size {n} and epsilon={EPSILON}")


sequences = ["".join(x) for x in product(("0","1"), repeat=n)]
probabilities = {seq:probability(seq) for seq in sequences}
print(f"{len(sequences)} sequences with total probability {sum(probabilities.values())}")

typical_set_min = 2**(-n*(binary_source_entropy(p)+EPSILON))
typical_set_max = 2**(-n*(binary_source_entropy(p)-EPSILON))
print(f"typical set limits: [{typical_set_min}, {typical_set_max}]")

typical_set = {seq:prob for seq, prob in probabilities.items() if typical_set_min<=prob<=typical_set_max}
typical_set_min_length = (1-EPSILON) * 2**(n*(binary_source_entropy(p)-EPSILON))
typical_set_max_length = 2**(n*(binary_source_entropy(p)+EPSILON))
typical_set_probability = sum(typical_set.values())
print(f"{len(typical_set)} typical sequences, expected between {typical_set_min_length} and {typical_set_max_length}")
print(f"total probability of the typical set {typical_set_probability}")

most_prob_seq = (None, 0)
for seq, prob in probabilities.items():
    if prob >= most_prob_seq[1]:
        most_prob_seq = (seq, prob)

print(f"the most probable sequence is {most_prob_seq[0]} with probability {most_prob_seq[1]} ({ones(most_prob_seq[0])} ones)")
if most_prob_seq in typical_set.items():
    print(f"the most probable sequence is in the typical set")
else:
    print(f"the most probable sequence is NOT in the typical set")


smallest_most_probable_set = {}
sorted_probs = dict(sorted(probabilities.items(), key=lambda item: item[1], reverse=True))
_sum = 0
for seq, prob in sorted_probs.items():
    if _sum < typical_set_probability:
        smallest_most_probable_set.update({seq: prob})
        _sum += prob
print(f"the size of the smallest set with probability {typical_set_probability} is {len(smallest_most_probable_set)}")

matches = 0
for seq in smallest_most_probable_set.keys():
    if seq in typical_set.keys():
        matches += 1
print(f"the typical set and the smallest equiprobable set have {matches} elements in common")
if most_prob_seq in smallest_most_probable_set.items():
    print(f"the most probable sequence is in the smallest most probable set")
else:
    print(f"the most probable sequence is NOT in the smallest most probable set >>>> SOMETHING IS WRONG <<<<")
