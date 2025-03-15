"""
Uniformly pick a probability p in [0, 1] and flip 2 coins with it. Measure how many times
both coins come up on the same side and how many times they differ.

The expected result is that they are the same 2/3 of the time.
"""

import random

num_of_experiments = 100_000_000 # increasing this will make the numbers converge more

same = 0
diff = 0

for i in range(num_of_experiments):
    p = random.random()
    trial_1 = 1 if random.random()<=p else 0
    trial_2 = 1 if random.random()<=p else 0
    
    if trial_1 == trial_2:
        same += 1
    else:
        diff += 1

    print(f"\r({i+1} trials): same: {same/(i+1)}, diff: {diff/(i+1)}", end="") # start with \r and end without \n to overwrite line
