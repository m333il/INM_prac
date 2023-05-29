from mpi4py import MPI
import numpy as np
import sys
import time

def norm(vec):
    local_norm = np.dot(vec, vec)
    global_norm = MPI.COMM_WORLD.allreduce(local_norm, MPI.SUM)
    return np.sqrt(global_norm)

def vec_sub(v1, v2, alpha):
    v1 -= alpha * v2

def dot(v1, v2):
    local_dot = np.dot(v1, v2)
    global_dot = MPI.COMM_WORLD.allreduce(local_dot, MPI.SUM)
    return global_dot

class Matrix:
    def __init__(self, row, col, my_rank, total_processes):
        self.row = row
        self.col = col
        self.my_rank = my_rank
        self.total_processes = total_processes
        self.A = np.random.rand(row, col)

        self.extra = 1
        if my_rank < col % total_processes:
            self.extra = 0
        
        self.offset = my_rank * row
        if self.extra:
            self.offset += col % total_processes
        
        self.dest = (my_rank + total_processes - 1) % total_processes
        self.source = (my_rank + 1) % total_processes

        self.v = np.zeros(row + self.extra)

        self.blocks_columns = [0 for _ in range(total_processes)]
        for i in range(total_processes):
            self.blocks_columns[i] = col // total_processes
            if (i + my_rank) % total_processes < col % total_processes:
                self.blocks_columns[i] += 1

    def matvec(self, x, res):
        np.copyto(self.v, x)
        res[:] = 0.0

        for k in range(self.total_processes):
            request_recieve = MPI.COMM_WORLD.Irecv(x, self.source, 21)
            request_send = MPI.COMM_WORLD.Isend(self.v, self.dest, 21)

            res[ : self.row] += self.A[:, self.offset : self.offset + self.blocks_columns[k]] @ self.v[ : self.blocks_columns[k]]
            self.offset = (self.offset + self.blocks_columns[k]) % self.col

            request_recieve.wait()
            request_send.wait()
            np.copyto(self.v, x)

    def save(self):
        np.savetxt(str(self.my_rank) + '.txt', self.A)
        

def arnoldi_iteration(A, k):
    n = A.row
    col = A.col
    my_rank = A.my_rank
    extra = A.extra

    v = np.zeros(n + extra)

    if my_rank == 0:
        h = np.zeros((k + 1) * k)
    else:
        h = np.array([])

    Q = np.zeros((k + 1, n + extra))
    np.copyto(Q[0], np.random.rand(n + extra))
    if extra:
        Q[0][-1] = 0.0
    Q0norm = norm(Q[0])
    Q[0] /= Q0norm
    
    for i in range(1, k + 1):
        A.matvec(Q[i - 1], v)

        for j in range(i):
            Qjv = dot(Q[j], v)
            if my_rank == 0:
                h[j * k + i - 1] = Qjv
            v -= Qjv * Q[j]
        
        vnorm = norm(v)
        if my_rank == 0:
            h[i * k + i - 1] = vnorm
        Q[i] = v / vnorm
    
    return h

def main():
    my_rank = MPI.COMM_WORLD.Get_rank()
    total_processes = MPI.COMM_WORLD.Get_size()

    n = int(sys.argv[1])
    k = int(sys.argv[2])

    my_amount = n // total_processes
    if my_rank < n % total_processes:
        my_amount += 1

    A = Matrix(my_amount, n, my_rank, total_processes)

    A.save()

    start = time.time()
    h = arnoldi_iteration(A, k)
    finish = time.time()

    if my_rank == 0:
        np.savetxt('h.txt', h)
        print(f'Amount of processes = {total_processes}, time = {finish - start}')

    MPI.Finalize()

if __name__ == '__main__':
    main()