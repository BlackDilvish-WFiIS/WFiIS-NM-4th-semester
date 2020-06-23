import matplotlib.pyplot as plt

X, Y, Z = [], [], []

for line in open(input('podaj plik: '), 'r'):
    values = [float(s) for s in line.split()]
    X.append(values[0])
    Y.append(values[1])
    #Z.append(values[2])

plt.plot(X, Y)
#plt.plot(X, Z)
#plt.yscale('log')
plt.show()