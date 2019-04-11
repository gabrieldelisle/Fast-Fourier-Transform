import matplotlib.pyplot as plt
import numpy as np



with open("sol.txt", 'r') as file :
	data = file.read()

data = [u.split(';') for u in data.split('\n')[:-1]]
t = [float(u[0]) for u in data]
x = [float(u[1]) for u in data]
X = [float(u[2]) for u in data]

A = t[0]
B = t[-1]
N = len(t)

plt.plot(t,x)
plt.show()


freqs = [k/(B-A) for k in range(N)]
plt.plot(freqs,X)
plt.show()


sp = np.fft.fft(x)
freqs = [k/(B-A) for k in range(N)]
plt.plot(freqs, np.abs(sp))
plt.show()