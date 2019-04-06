import matplotlib.pyplot as plt

with open("sol.txt", 'r') as file :
	data = file.read()

data = [u.split(';') for u in data.split('\n')[:-1]]
t = [float(u[0]) for u in data]
x = [float(u[1]) for u in data]
X = [float(u[2]) for u in data]

plt.plot(t,x)
plt.show()
plt.plot(t,X)
plt.show()