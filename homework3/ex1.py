import random

rates = [1, 1/2, 1/3, 1/4, 1/5, 1/8]
p = 1/10

def random_N_R_code_distinct(n, r):
    """
    Create a random binary code with codeword length N and rate R, with
    distinct codewords.
    """
    num_of_codewords = 2**int(n*r)

    # random.sample() fails for possible codewords = 2**64 with OverflowError: Python int too large to convert to C ssize_t
    randomly_picked_codewords = set()
    while len(randomly_picked_codewords) < num_of_codewords:
        randomly_picked_codewords.add(random.randint(0, 2**n))

    return [[1 if bit=='1' else 0 for bit in f"{codeword:0{n}b}"] for codeword in randomly_picked_codewords]


def random_N_R_code(n, r):
    """
    Create a random binary code with codeword length N and rate R.
    """
    # we code NR bits
    num_of_codewords = 2**int(n*r)

    # we create 2**NR codewords at random
    codewords = [[0 for _ in range(n)] for _ in range(num_of_codewords)]
    for i in range(num_of_codewords):
        for j in range(n):
            x_i = random.randint(0, 1)
            codewords[i][j] = x_i

    return codewords

def channel(x, p):
    output = [0] * len(x)
    for i, b in enumerate(x):
        output[i] = b if random.random() >= p else int(not b)
    return output

def ML_decode(received, inputs):
    min_dist = sum([1 for x, y in zip(received, inputs[0]) if x!=y])
    estimate = inputs[0]
    for input in inputs[1:]:
        hamming_dist = sum([1 for x, y in zip(received, input) if x!=y])
        if hamming_dist < min_dist:     # assuming no equalities
            min_dist = hamming_dist
            estimate = input
    return estimate

for rate in rates:

    unit = int(1/rate)

    bit_errors = [0 for _ in range(10)]

    for n in range(unit, 10*unit+1, unit):
        codewords = random_N_R_code_distinct(n, rate)
        print(f"rate {rate}, n {n}, unit {unit}, {len(codewords)} ({2**int(n*rate)}) codewords, {len(set([''.join(map(str, word)) for word in codewords]))} distinct")

        for i in range(100_000):
            print(f"experiment {i}", end="\r", flush=True)

            cw_idx = random.randint(0, len(codewords)-1)
            cw = codewords[cw_idx]

            channeled = channel(cw, p)
            decoded = ML_decode(channeled, codewords)
            _bit_errors = sum([1 for x, y in zip(cw, decoded) if x!=y])

            bit_errors[int(n/unit)-1] += _bit_errors
    print(f"bit errors for rate {rate}: {bit_errors}")
