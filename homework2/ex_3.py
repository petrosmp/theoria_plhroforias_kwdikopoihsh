import heapq
from itertools import product
from math import log2

def binary_source_entropy(p):
    return -p*log2(p) - (1-p)*log2(1-p) 

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

    return codebook, [len(x) for x in codebook]

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
    
    return shannon_fano_codebook(probs)

def huffman_codebook(probs):
    heap = [[prob, [i, ""]] for i, prob in enumerate(probs)]
    heapq.heapify(heap)


    while len(heap) > 1:
        a = heapq.heappop(heap)
        b = heapq.heappop(heap)

        for pair in a[1:]:
            pair[1] = '0' + pair[1]
        for pair in b[1:]:
            pair[1] = '1' + pair[1]
 
        heapq.heappush(heap, [a[0] + b[0]] + a[1:] + b[1:])

    # only 1 element remains at the end
    codes = heap[0][1:]

    codebook = [None for _ in probs]
    for symbol, code in codes:
        codebook[symbol] = code

    return codebook, [len(x) for x in codebook]

def huffman_codebook_seq(n, p):
    """
    Generate a Huffman codebook for all binary sequences of length n, where 1 appears with probability p.
    """
    probs = []
    sequences = ["".join(x) for x in product(("0","1"), repeat=n)]
    for seq in sequences:
        ones = sum([1 for i in seq if i=="1"])
        probs.append(p**ones * (1-p)**(n-ones))

    return huffman_codebook(probs)

def encoding(X, n, codebook, code_lengths):
    return codebook[X]

def decoding(Y, p, n, codebook, code_lengths):
    for i in range(2**n):
        if codebook[i] == Y:
            return i


hcodebook, h_lens = huffman_codebook_seq(4, 0.1)
sfcodebook, sf_lens = shannon_fano_codebook_seq(4, 0.1)

print(f"huffman codebook:\n{hcodebook}\n\n\n")
print(f"shannon fano codebook:\n{sfcodebook}")


raw = 0b0011
enc = encoding(raw, 4, hcodebook, h_lens)
dec = decoding(enc, 0.1, 4, hcodebook, h_lens)

print(raw)
print(enc)
print(dec)
