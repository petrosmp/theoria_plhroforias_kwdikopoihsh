"""
A Discrete Memoryless Source (DMS) produces a sequence {X_i}, i = 1,2,3,...
where P(X_i=1)=p for all i. Its output is being processed in blocks of n bits.

1) Where are typical sequences observed
"""
import matplotlib.pyplot as plt
from math import log2

EPSILON = 0.05

p = 0.45
n = 15

def prob_of_block_with_index(p:int, i: int):
    """
    Return the probability of a block that can be converted to the decimal number i when
    the probability of a 1 is p.
    """

    # convert i to binary
    binary_repr = format(i, 'b')

    # count ones
    ones = sum([int(x) for x in binary_repr])

    # find probability
    return pow(p, ones)*pow(1-p,n-ones)

probabilities = []

for i in range(2**n):
    probabilities.append(prob_of_block_with_index(p, i))

def binary_source_entropy(p):
    return -p*log2(p) - (1-p)*log2(1-p) 

probabilities.sort(reverse=True)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('symlog')
ax.set_xscale('symlog')
ax.set_xlim(-0.25*10**4, 3.5*10**4)
line, = ax.plot(probabilities, color='blue', lw=2)
ax.axhline(pow(2, -n*binary_source_entropy(p)))
ax.axhline(pow(2, -n*(binary_source_entropy(p)-EPSILON)))
ax.axhline(pow(2, -n*(binary_source_entropy(p)+EPSILON)))
ax.set_xlim(-0.25*10**4, 3.5*10**4)
plt.show()
