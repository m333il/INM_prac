import numpy as np

A = np.loadtxt('A.txt')
h = np.loadtxt('h.txt')

eigA = np.linalg.eig(A)[0]
eigh = np.linalg.eig(h)[0]
print(eigA)
print(eigh)