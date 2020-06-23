import matplotlib.pyplot as plt

X, Y, Z = [], [], []

with open("data15.dat", "r") as f:
    for line in f:
        dat = line.split()
        X.append(float(dat[0]))
        Y.append(float(dat[1]))
        Z.append(float(dat[2]))

# plt.ylim(-3, 3)
plt.plot(X, Y)
plt.plot(X, Z)
plt.show()
