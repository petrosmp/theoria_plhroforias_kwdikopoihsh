from math import log2, comb
import matplotlib.pyplot as plt

EPSILON = 0.05

def binary_source_entropy(p):
    return -p*log2(p) - (1-p)*log2(1-p) 


p = 0.45
H  = binary_source_entropy(p)

A_ns = []
upper_bounds = []
lower_bounds = []

for n in range(10, 1000+1, 10):
    
    # there are n choose np typical sequences
    A_n = comb(n, int(n*p))

    upper_bound = pow(2, n*(H+EPSILON))
    lower_bound = (1-EPSILON) * pow(2, n*(H-EPSILON))

    A_ns.append(A_n)
    upper_bounds.append(upper_bound)
    lower_bounds.append(lower_bound)

    print(f"{n}")



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('symlog')
line, = ax.plot(n, A_ns, color='blue')
line, = ax.plot(n, upper_bounds, color='red')
line, = ax.plot(n, lower_bounds, color='green')

plt.show()