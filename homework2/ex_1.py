"""
Source with alphabet X={1,2,3,4} with probabilities 1/2, 1/4, 1/8, 1/8.

Huffman code:
    1 -> 0
    2 -> 10
    3 -> 110
    4 -> 111
"""

import random
from collections import Counter
import matplotlib.pyplot as plt

X = [1,2,3,4]
entropy = 1.75
n = 100
codebook = ["0", "10", "110", "111"]
num_of_experiments = 1000

def huffman(raw, codebook):
    encoded = ""
    for elem in raw:
        encoded += codebook[elem-1]
    return encoded

def dehuffman(encoded, codebook):
    decoded = []
    tmp = ""
    for char in encoded:
        tmp += char
        if tmp in codebook:
            decoded.append(codebook.index(tmp)+1)
            tmp = ""
    return decoded

errors = 0
avg_code_word_lengths = []
d_5 = []
d_6 = []
joint = []
for _ in range(num_of_experiments):

    raw = random.choices(X, k=n, weights=[1/2, 1/4, 1/8, 1/8])
    encoded = huffman(raw, codebook)
    decoded = dehuffman(encoded, codebook)

    avg_code_word_length = len(encoded) / n
    avg_code_word_lengths.append(avg_code_word_length)

    print(encoded)
    print(decoded)

    d_5.append(encoded[4])
    d_6.append(encoded[5])
    joint.append(encoded[4]+encoded[5])

    for i, (r, d) in enumerate(zip(raw, decoded)):
        if r != d:
            print(f"error at position {i}, raw [{r}] != decoded [{d}]")
            errors += 1

print(f"errors: {errors}")


avg_code_word_lengths = Counter(avg_code_word_lengths)
avg_code_word_lengths = dict(sorted(avg_code_word_lengths.items(), key=lambda item: item[0], reverse=True))

d_5 = dict(sorted(Counter(d_5).items(), key=lambda item: item[0]))
d_6 = dict(sorted(Counter(d_6).items(), key=lambda item: item[0]))
joint = dict(sorted(Counter(joint).items(), key=lambda item: item[0]))
print(d_5)
print(d_6)
print(joint)



plt.plot(avg_code_word_lengths.keys(), avg_code_word_lengths.values())
plt.title("Distribution of average code word length")
plt.axvline(x=entropy, color="black", linestyle="dashed", label="entropy")
plt.legend()

plt.figure()
plt.subplot(2, 2, 1)
plt.title("Distribution of d_5")
plt.bar(d_5.keys(), d_5.values())
plt.subplot(2, 2, 2)
plt.title("Distribution of d_6")
plt.bar(d_6.keys(), d_6.values(), color="orange")



plt.subplot(2, 1, 2)
plt.title("Joint distribution of d_5 and d_6")
plt.bar(joint.keys(), joint.values(), color="green")


plt.show()
