import numpy as np

A = np.loadtxt('A.txt')
h = np.loadtxt('h.txt')

eigA = np.linalg.eig(A)[0]
eigA = eigA[np.argsort(np.abs(eigA))]
print(eigA)
eigh = np.linalg.eig(h)[0]
eigh = eigh[np.argsort(np.abs(eigh))]
print(eigh)