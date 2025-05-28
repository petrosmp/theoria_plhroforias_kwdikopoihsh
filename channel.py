from math import log2, isclose

def H(p):
    return -p*log2(p) - (1-p)*log2(1-p) 

x_distr = [3/10, 7/10]
H_x = H(x_distr[0])
print(f"x distribution {x_distr} (sum {sum(x_distr)})")
print(f"H(X) = {H_x}")

channel_transition_matrix = [
    [2/3, 1/3],
    [1/3, 2/3]
]
print(f"channel: {channel_transition_matrix} (sum {sum([sum(x) for x in channel_transition_matrix])})")

y_distr = [0, 0]
for i in range(len(x_distr)):
    for j in range(len(y_distr)):
        y_distr[j] += x_distr[i] * channel_transition_matrix[i][j] # prob to be in x[i] and end up in y[j]
H_y = H(y_distr[0])
print(f"y distribution {y_distr} (sum {sum(y_distr)})")
print(f"H(Y) = {H_y}")

joint_distr = []
for i in range(len(x_distr)):
    tmp = []
    for j in range(len(y_distr)):
        tmp.append(x_distr[i]*channel_transition_matrix[i][j])
    joint_distr.append(tmp)
print(f"joint distribution: {joint_distr} (sum {sum([sum(x) for x in joint_distr])})")

mutual_information = 0
for i in range(len(x_distr)):
    for j in range(len(y_distr)):
        p_xy = joint_distr[i][j]
        p_x = x_distr[i]
        p_y = y_distr[j]
        mutual_information += p_xy * log2(p_xy/(p_x*p_y))

print(f"mutual information: {mutual_information}")

reverse_channel = []
for i in range(len(x_distr)):
    tmp = []
    for j in range(len(y_distr)):
        tmp.append(joint_distr[i][j] / y_distr[j])
    reverse_channel.append(tmp)
print(f"reverse channel: {reverse_channel} (sum {sum([sum(x) for x in reverse_channel])})")


H_y_given_x = 0
for i in range(len(x_distr)):
    for j in range(len(y_distr)):
        H_y_given_x -= joint_distr[i][j] * log2(channel_transition_matrix[i][j])
print(f"H(Y|X) = {H_y_given_x}")
print(f"H(Y) - H(Y|X) = {H_y - H_y_given_x} ({'not ' if not isclose(H_y-H_y_given_x, mutual_information) else ''}equal to mutual information)")

H_x_given_y = 0
for i in range(len(x_distr)):
    for j in range(len(y_distr)):
        H_x_given_y -= joint_distr[i][j] * log2(reverse_channel[i][j])
print(f"H(X|Y) = {H_x_given_y}")
print(f"H(X) - H(X|Y) = {H_x - H_x_given_y} ({'not ' if not isclose(H_x-H_x_given_y, mutual_information) else ''}equal to mutual information)")

