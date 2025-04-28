from math import ceil, log2
from itertools import product

def create_tree(max_depth):
    """
    The "tree" is an array of nodes represented by pairs (left, right) specifying the index
    of its left and right child respectively.

    NULL is represented by -1.
    """
    tree = []
    
    # create internal levels
    for i in range(1, 2**(max_depth)-1, 2):
        tree += [[i, i+1]]

    # create bottom level
    for i in range(0, 2**(max_depth)-2**(max_depth-1)):
        tree += [[-1, -1]]

    return tree

def shannon_fano(probs):

    lengths = []
    for prob in probs:
        lengths.append(ceil(-log2(prob)))

    lengths = sorted(lengths)

    max_length = max(lengths)
    tree = create_tree(max_length+1)

    codes = []
    for length in lengths:
        code = []
        cur = 0
        step = 0
        while step < length:
            if tree[cur][0] != -1:
                if step == length-1:
                    tree[cur][0] = -1
                cur = tree[cur][0]
                step += 1
                code += [0]
            elif tree[cur][1] != -1:
                if step == length-1:
                    tree[cur][1] = -1
                cur = tree[cur][1]
                step += 1
                code += [1]
            else:
                parent = ceil((cur-2)/2)
                if parent*2 + 1 == cur: # left
                    tree[parent][0] = -1
                else:   # right
                    tree[parent][1] = -1
                cur = 0
                step = 0
                code = []
        codes.append(code)
    
    return codes

def shannon_fano_generation(n, p):
    probs = []
    sequences = ["".join(x) for x in product(("0","1"), repeat=n)]
    for seq in sequences:
        ones = sum([1 for i in seq if i=="1"])
        probs.append(p**ones * (1-p)**(n-ones))

    return shannon_fano(probs)

n = 10
p = 0.1
codes = shannon_fano_generation(n ,p)
for i, code in enumerate(codes):
    print(f"code {i}: {code} (len {len(code)})")
    
codes = shannon_fano([3/8, 3/8, 1/8, 1/8])
print(codes)
codes = shannon_fano([1/2, 1/4, 1/8, 1/8])
print(codes)
