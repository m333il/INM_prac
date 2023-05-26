import numpy as np

amount = 4
n = 1000

As = []

for i in range(amount):
    As.append(np.loadtxt(str(i) + '.txt'))

A = np.array(As[0])

for i in range(1, amount):
    A = np.concatenate((A, As[i]))
A = A.reshape(n, n)
#print(A)
#A = np.array(As)
#A = A.reshape(10, 10);

eigA = np.linalg.eig(A)[0]
eigA = np.sort_complex(eigA)
#print(eigA)

H = np.loadtxt('h.txt').reshape(n, n)

eigH = np.linalg.eig(H)[0]
eigH = np.sort_complex(eigH)
#print(eigH)

wr = np.loadtxt('wr.txt')
wi = np.loadtxt('wi.txt')

eigh = wr + wi * 1j

eigh = np.sort_complex(eigh)
#print(eigh)

print(np.abs(eigA - eigh))

print(f'error = {np.linalg.norm(eigA - eigh)}')
print(f'error = {np.linalg.norm(eigH - eigh)}')
