import matplotlib.pyplot as plt
import numpy as np

# read data from output file
with open("sol.txt", "r") as file:
    data = file.read()

data = [u.split(";") for u in data.split("\n")[:-1]]
t = [float(u[0]) for u in data]
x = [float(u[1]) for u in data]
X = [float(u[2]) for u in data]

A = t[0]
B = t[-1]
N = len(t)

# display initial array
plt.plot(t, x)
plt.show()

# compare our result with numpy
freqs = [k / (B - A) for k in range(N)]
sp = np.fft.fft(x)
sp = np.abs(sp)
plt.plot(freqs, X, label="ours")
plt.plot(freqs, sp, label="numpy")
plt.legend()
plt.show()
print(np.allclose(X, sp))
