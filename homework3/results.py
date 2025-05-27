import matplotlib.pyplot as plt

bit_errors = [
    [9961, 22140, 30998, 40865, 51058, 60237, 70150, 80492, 89672, 99759],
    [20014, 31194, 40091, 60474, 70752, 86504, 98936, 110828, 123494, 136683],
    [19856, 26744, 39962, 44417, 46274, 46188, 49450, 51736, 48632, 49090],
    [19938, 23969, 22111, 24253, 23342, 20072, 17544, 14577, 12626, 10858],
    [8571, 23226, 13574, 12833, 8845, 7303, 5012, 3592, 2613, 1392],
    [8448, 6909, 2146, 1083, 888, 193, 0, 0, 28, 0],
]

rates = [1, 1/2, 1/3, 1/4, 1/5, 1/8]
p = 1/10

ber = []

for i in range(len(rates)):
    _ber = []
    errors = bit_errors[i]
    rate = rates[i]

    unit = int(1/rate)

    lengths = list(range(unit, 10*unit+1, unit))

    for j, length in enumerate(lengths):
        bit_error_rate = errors[j]/(length*100_000)
        _ber.append(bit_error_rate)
    ber.append(_ber)

for i in range(len(ber)):
    plt.plot(list(range(1, 11)), ber[i], label=f"R = {round(rates[i], 3)}")


plt.title("Evolution of Bit Error Rate with Codeword Length")
plt.xlabel("Sequence Length (in bits, x$R^{-1}$)")
plt.ylabel("Bit Error Rate")
plt.legend()

plt.show()
