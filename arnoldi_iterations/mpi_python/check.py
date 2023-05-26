import numpy as np
import sys

amount = sys.argv[1]
n = sys.argv[2]

As = []

for i in range(amount):
    As.append(np.loadtxt(str(i) + '.txt'))

A = np.array(As[0])

for i in range(1, amount):
    A = np.concatenate((A, As[i]))
A = A.reshape(n, n)

eigA = np.linalg.eig(A)[0]
eigA = np.sort_complex(eigA)

H = np.loadtxt('h.txt')[: n * n].reshape(n, n)

eigH = np.linalg.eig(H)[0]
eigH = np.sort_complex(eigH)

print(f'error = {np.linalg.norm(eigA - eigH)}')