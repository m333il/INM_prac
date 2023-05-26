#include "mpi.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <algorithm>
#include <chrono>
#include <lapacke.h>
#include <cstring>
#include <random>

template<class T>
T norm(T *vec, int n) {
    double local = 0.0;
    for (unsigned i = 0; i < n; ++i)
        local += vec[i] * vec[i];
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(global);
}

template<class T>
void normalize(T *vec, int n) {
    double vnorm = norm(vec, n);
    for (int i = 0; i < n; ++i) {
        vec[i] /= vnorm;
    }
}

template<class T>
void vec_sub(T *v1, T *v2, T alpha, int n) {
    for (int i = 0; i < n; ++i) {
        v1[i] -= alpha * v2[i];
    }
}

template<class T>
T dot(T *v1, T *v2, int n) {
    double local = 0.0;
    for (int i = 0; i < n; ++i) {
        local += v1[i] * v2[i];
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

class Matrix
{
public:
    int col, row;
    int my_rank, total_processes;
    int extra, offset;
    int dest, source;
    double *A;
    double *v;
    MPI_Request requests[2];
    MPI_Status statuses[2];

    Matrix(int row, int col, int my_rank, int total_processes) : 
        col(col), row(row), my_rank(my_rank), total_processes(total_processes)
    {
        A = new double[row * col];
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution <>floatDist(-10.0, 10.0);
        for (int i = 0; i < row * col; ++i) {
            A[i] = floatDist(gen);
        }

        // we need a space for extra element in processes that don't have one
        // because we need to send/recieve arrays of same length
        extra = my_rank < col % total_processes ? 0 : 1;
        // offset of local A in whole A
        offset = my_rank * row;
        if (extra)
            offset += col % total_processes;

        // next process
        dest = (my_rank + total_processes - 1) % total_processes;

        // prev process
        source = (my_rank + 1) % total_processes;

        // array to store current part of vector
        v = new double[row + extra];
    }

    void matvec(double *x, double *res) {
        memcpy(v, x, sizeof(double) * (row + extra));
        std::fill(res, res + row, 0.0);

        for (int k = 0; k < total_processes; ++k) {
            int col_block = col / total_processes;
            if ((k + my_rank) % total_processes < col % total_processes) {
                ++col_block;
            }
            // send our data to next process and recieve new data from previous process
            MPI_Irecv(x, row + extra, MPI_DOUBLE, source, 21, MPI_COMM_WORLD, &requests[0]);
            MPI_Isend(v, row + extra, MPI_DOUBLE, dest, 21, MPI_COMM_WORLD, &requests[1]);
            // while getting new data - perform work on current data
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < col_block; ++j) {
                    res[i] += A[i * col + (j + offset)] * v[j];
                }
            }
            offset = (offset + col_block) % col;
            // wait for getting new data
            MPI_Waitall(2, requests, statuses);
            memcpy(v, x, sizeof(double) * (row + extra));
        }
    }

    uint32_t rows() const {
        return row;
    }

    uint32_t columns() const {
        return col;
    }

    uint32_t rank() const {
        return my_rank;
    }

    uint32_t get_extra() const {
        return extra;
    }

    void save() const {
        std::ofstream fout(std::to_string(my_rank) + ".txt");
        for (int i = 0; i < row * col; ++i) {
            fout << A[i] << " ";
        }
        fout.close();
    }

    ~Matrix() {
        delete[] A;
        delete[] v;
    }
};

template<class T, class Matrix>
T *arnoldi_iteration(Matrix &A, int k) {
    uint32_t n = A.rows();
    uint32_t col = A.columns();
    uint32_t my_rank = A.rank();

    int extra = A.get_extra();
    T *Q = new T[(n + extra) * (k + 1)];
    std::fill(Q, Q + (n + extra) * (k + 1), 0.0);
    T *h = nullptr;

    if (my_rank == 0) {
        // h is stored in only one process
        // because it doesn't take part in calculations as a vector
        h = new T[(k + 1) * k];
        std::fill(h, h + (k + 1) * k, 0.0);
    }

    T *v = new T[n];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution <>floatDist(-10.0, 10.0);
    for (int i = 0; i < n; ++i) {
        Q[i] = floatDist(gen);
    }
    normalize(Q, n);

    double vnorm;
    for (int i = 1; i <= k; ++i) {
        A.matvec(Q + (i - 1) * (n + extra), v);

        for (int j = 0; j < i; ++j) {
            double Qjv = dot(Q + j * (n + extra), v, n);
            if (my_rank == 0) {
                h[j * k + i - 1] = Qjv;
            }
            vec_sub(v, Q + j * (n + extra), Qjv, n);
        }
        vnorm = norm(v, n);
        if (my_rank == 0) {
            h[i * k + i - 1] = vnorm;
        }
        for (int j = 0; j < n; ++j) {
            Q[i * (n + extra) + j] = v[j] / vnorm;
        }
    }

    delete[] v;
    delete[] Q;

    return h;
}

int main(int argc, char *argv[]) {    
    MPI_Init(&argc, &argv);

    int my_rank, total_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_processes);

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);

    int my_amount = n / total_processes;
    if (my_rank < n % total_processes) {
        ++my_amount;
    }

    Matrix A(my_amount, n, my_rank, total_processes);
    A.save();

    double *h = arnoldi_iteration<double, Matrix>(A, k);
    if (my_rank == 0) {
        std::ofstream fouth("h.txt");
        for (int i = 0; i < k * k; ++i) {
            fouth << h[i] << " ";
        }
        fouth.close();
        double *wr = new double[k];
        double *wi = new double[k];
        double *Z = new double[k * k];
        int info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'E', 'N', k, 1, k, h, k, wr, wi, Z, k);
        std::ofstream foutwi("wi.txt");
        std::ofstream foutwr("wr.txt");
        
        for (int i = 0; i < k; ++i) {
            foutwi.precision(17);
            foutwr.precision(17);
            foutwi << std::fixed << wi[i] << " ";
            foutwr << std::fixed << wr[i] << " ";
        }
        foutwi.close();
        foutwr.close();
        delete[] Z;
        delete[] wr;
        delete[] wi;
    }
    if (my_rank == 0) {
        delete[] h;
    }
    MPI_Finalize();
}
