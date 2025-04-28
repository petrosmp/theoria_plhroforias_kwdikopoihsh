import random
import heapq
from collections import Counter
import matplotlib.pyplot as plt

transition_matrix = [
   [1/2,  1/8,  1/8,  1/8,  1/8],
   [1/4,  1/8,  1/16, 1/16, 1/2],
   [1/4,  1/8,  1/8,  1/4,  1/4],
   [1/8,  0,    1/2,  1/4,  1/8],
   [0,    1/2,  1/4,  1/4,  0],
]

stationary_distr = [1/5, 1/5, 1/5, 1/5, 1/5]

num_of_experiments = 1000
n = 100

def huffman_codebook(probs):
    heap = [[prob, [i+1, ""]] for i, prob in enumerate(probs)]
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
        codebook[symbol-1] = code

    return codebook

def generate_sequence(length, transition_matrix, initial_distr):
    
    state = random.choices(range(len(transition_matrix)), initial_distr)[0]
    sequence = [state+1]

    for _ in range(length-1):
        state = random.choices(range(len(transition_matrix)), transition_matrix[state])[0]
        sequence.append(state+1)

    return sequence

def huffman(raw, codebook):
    encoded = ""
    for elem in raw:
        encoded += codebook[elem-1]
    return encoded

def dehuffman(encoded, codebook):
    decoded = []
    tmp = ""
    for char in encoded:
        tmp += str(char)
        if tmp in codebook:
            decoded.append(codebook.index(tmp)+1)
            tmp = ""
    return decoded

def markov_huffman(raw, transition_matrix, initial_prob):
    encoded = ""
    codebook = huffman_codebook(initial_prob)
    for state in raw:
        encoded += codebook[state-1]
        codebook = huffman_codebook(transition_matrix[state-1])
    return encoded

def markov_dehuffman(encoded, transition_matrix, initial_prob):
    decoded = []
    codebook = huffman_codebook(initial_prob)
    tmp = ""
    for char in encoded:
        tmp += char
        if tmp in codebook:
            symbol = codebook.index(tmp)+1
            decoded.append(symbol)
            tmp = ""
            codebook = huffman_codebook(transition_matrix[symbol-1])
    return decoded

initial_prob = [random.random() for _ in range(len(transition_matrix))]
initial_prob = [x/sum(initial_prob) for x in initial_prob]


stationary_codebook = huffman_codebook(stationary_distr)

stationary_avg_code_word_lengths = []
conditional_avg_code_word_lengths = []

stationary_entropy = 2.32
entropy_rate = 1.875 
stationary_errors = 0
conditional_errors = 0

for _ in range(num_of_experiments):

    initial_prob = [random.random() for _ in range(len(transition_matrix))]
    initial_prob = [x/sum(initial_prob) for x in initial_prob]

    seq = generate_sequence(n, transition_matrix, initial_prob)

    stationary_enc = huffman(seq, stationary_codebook)
    conditional_enc = markov_huffman(seq, transition_matrix, initial_prob)

    stationary_dec = dehuffman(stationary_enc, stationary_codebook)
    conditional_dec = markov_dehuffman(conditional_enc, transition_matrix, initial_prob)

    for r, s, c in zip(seq, stationary_dec, conditional_dec):
        if s != r:
            stationary_errors += 1
        if c != r:
            conditional_errors += 1

    stationary_avg_code_word_lengths.append(len(stationary_enc)/n)
    conditional_avg_code_word_lengths.append(len(conditional_enc)/n)

print(f"errors using stationary huffman: {stationary_errors}")
print(f"errors using conditional huffman: {conditional_errors}")

stationary_avg_code_word_lengths = Counter(stationary_avg_code_word_lengths)
conditional_avg_code_word_lengths = Counter(conditional_avg_code_word_lengths)
stationary_avg_code_word_lengths = dict(sorted(stationary_avg_code_word_lengths.items(), key=lambda item: item[0], reverse=True))
conditional_avg_code_word_lengths = dict(sorted(conditional_avg_code_word_lengths.items(), key=lambda item: item[0], reverse=True))

plt.plot(stationary_avg_code_word_lengths.keys(), stationary_avg_code_word_lengths.values(), label="avg stationary length")
plt.plot(conditional_avg_code_word_lengths.keys(), conditional_avg_code_word_lengths.values(), label="avg conditional length")
plt.title("Distribution of average code word length")
plt.axvline(x=stationary_entropy, color="black", linestyle="dashed", label="stationary entropy")
plt.axvline(x=entropy_rate, color="red", linestyle="dashed", label="entropy rate")
plt.legend()

plt.show()
