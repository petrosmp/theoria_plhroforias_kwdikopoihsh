"""
A Discrete Memoryless Source (DMS) produces a sequence {X_i}, i = 1,2,3,...
where P(X_i=1)=p for all i. Its output is being processed in blocks of n bits.

This script plots the probabilities of blocks w.r.t. their indices, in decreasing
order.


We wish to plot the probabilities of blocks with respect to their indices,
sorted from largest to smallest (the index of a block is its decimal
representation).
The fact that we are sorting them effectively means that the resulting
plot will be a piecewise function looking like a staircase, where the
values are the probabilities of a block with a number of 1's and the
length of each "step" will be the number of ways in which it is possible
to have that number of 1's in a block.
In other words, for i ranging from 0 to n (included), we want to plot
each (p^i)*(1-p)^(n-i), (n choose i) times. An effective way to do this
without calculating probabilities for all 2^n possible indices is to
create an array with (n choose i) elements of each probability. For large
values of n (i.e. 50) however this array gets prohibitively big (2^50 is
close to the number of grains of sand on earth, we cannot possibly fit
that many values in our RAM), and thus we employ a clever sampling technique:
  - we calculate the range which each probability covers. for example,
    the first (n choose 1) numbers all correspond to the probability of
    a single one in a string, the next (n choose 2) correspond to the
    probability of 2 ones in a string etc
  - we divide [0, 2^n] in powers of 10
  - for each power of 10, we take 100 equally spaced samples
  - we match each sample to probability of the range it belongs to
The cheatsheet is of the form [(cumulative n choose i, probability of i ones
in a block)] so that by searching it in reverse we can find the range in
which a given number falls and thus match it with a probability in O(n) time.

Using this technique we effectively reduce the number of samples we
plot from 2^n to O(n), allowing us to handle otherwise unapproachable
values of n.
"""


import matplotlib.pyplot as plt
from math import log2, log10, comb
import time


EPSILON = 0.05


p = 0.45
n = 50
samples_per_power_of_ten = 10000
xstart = int(0)
xend = int(12e14)

def binary_source_entropy(p):
    return -p*log2(p) - (1-p)*log2(1-p) 

start = time.time()

cheatsheet = []
counter = 0
for i in range(n+1):
    n_choose_i = comb(n, i)
    prob_i_ones = pow(p, i)*pow(1-p, n-i)
    counter += n_choose_i
    cheatsheet.insert(0, (counter, prob_i_ones))

def efficient(x):
    for counter, prob in cheatsheet:
        if x >= counter:
            return prob
    return cheatsheet[-1][1]

start = xstart
x = []
for i in range(int(log10(xend))):
    power_of_10_increment = log10(samples_per_power_of_ten)
    until = min(int(10**(i+power_of_10_increment)), xend+1)
    # go from start to 10^i + (desired_samples)*(10^i) in increments of 10^i
    x.extend(list(range(start, until, int(10**(i)))))
    start = until

print(f"max is {max(efficient(i) for i in x)}")

end = time.time()
print(f"time: {end-start} s")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('symlog')
ax.set_xscale('symlog')
ax.set_xlim(xstart, xend)
line, = ax.plot(x, [efficient(i) for i in x], color='blue', lw=2)
ax.axhline(pow(2, -n*binary_source_entropy(p)))
ax.axhline(pow(2, -n*(binary_source_entropy(p)-EPSILON)))
ax.axhline(pow(2, -n*(binary_source_entropy(p)+EPSILON)))

plt.show()