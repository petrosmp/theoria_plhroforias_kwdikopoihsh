from math import log2, comb
import matplotlib.pyplot as plt

EPSILON = 0.05

def binary_source_entropy(p):
    return -p*log2(p) - (1-p)*log2(1-p) 


p = 0.45
H = binary_source_entropy(p)



probs = []
lower_bounds = []

for n in range(10, 1000+1, 10):
    
    prob = comb(n, int(n*p))*pow(p, n*p)*pow(1-p, n-n*p)

    lower_bound = 1-EPSILON

    print(f"prob {prob}")

    probs.append(prob)
    lower_bounds.append(lower_bound)



n = list(range(10, 1000+1, 10))
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(n, probs, color='blue')
line, = ax.plot(n, lower_bounds, color='green')

plt.show()