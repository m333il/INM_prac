#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <lapacke.h>

template<class T>
T trand(T a, T b) {
    double f = (double)rand() / RAND_MAX;
    return (T)(a + f * (b - a));
}

template<class T>
class Matrix
{
private:
    const uint32_t n;
    T *A;
    const T a, b;
public:
    Matrix(uint32_t size, T a = -10.0, T b = 10.0) : 
        n(size), a(a), b(b)
    {
        A = new T[n * n];
        for (uint32_t i = 0; i < n * n; ++i) {
            A[i] = trand(a, b);
        }
    }

    uint32_t order() const {
        return n;
    }

    void save(const std::string &filename) {
        std::ofstream A_file("A.txt");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A_file << A[i * n + j] << " ";
            }
            A_file << std::endl;
        }
        A_file.close();
    }

    void matvec(T *x, T *res) const {
        uint32_t i, j;
        for (i = 0; i < n; ++i) {
            res[i] = A[i * n] * x[0];
            for (j = 1; j < n; ++j) {
                res[i] += A[i * n + j] * x[j];
            }   
        }
    }

    ~Matrix() {
        delete[] A;
    }
};

template<class T>
T norm(T *vec, int n) {
    T res = 0.0;
    for (unsigned i = 0; i < n; ++i)
        res += vec[i] * vec[i];
    return sqrt(res);
}

template<class T>
void normalize(T *vec, int n) {
    T vnorm = norm(vec, n);
    for (int i = 0; i < n; ++i) {
        vec[i] /= vnorm;
    }
}

template<class T>
T dot(T *v1, T *v2, int n) {
    T res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

template<class T>
void vec_sub(T *v1, T *v2, T alpha, int n) {
    for (int i = 0; i < n; ++i) {
        v1[i] -= alpha * v2[i];
    }
}

template<class T>
void arnoldi_iteration(Matrix<T> &A, int k, T *Q, T *h) {
    //A - n x n array
    //b - initial vector, length n
    //k - amount of iterations

    const uint32_t n = A.order();
    double eps = 1e-12;
    T *v = new T[n];

    for (int i = 0; i < n; ++i) {
        Q[i] = trand(-10.0, 10.0);
    }
    normalize(Q, n);

    for (int i = 1; i <= k; ++i) {
        A.matvec(Q + (i - 1) * n, v);

        for (int j = 0; j < i; ++j) {
            h[j * k + i - 1] = dot(Q + j * n, v, n);
            vec_sub(v, Q + j * n, h[j * k + i - 1], n);
        }
        h[i * k + i - 1] = norm(v, n);
        if (h[i * k + i - 1] <= eps) {
            break;
        }
        for (int j = 0; j < n; ++j) {
            Q[i * n + j] = v[j] / h[i * k + i - 1];
        }
    }
    delete[] v;
}



int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    Matrix<double> A(n);

    double *Q = new double[(k + 1) * n];
    double *h = new double[(k + 1) * k];
    std::fill(Q, Q + (k + 1) * n, 0);
    std::fill(h, h + (k + 1) * k, 0);
    double *wr = new double[k];
    double *wi = new double[k];
    double *Z = new double[k * k];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    arnoldi_iteration(A, k, Q, h);
    
    int info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'E', 'N', k, 1, k, h, k, wr, wi, Z, k);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time for n = " << n << " is " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
    delete[] Z;
    delete[] Q;
    delete[] h;
    delete[] wr;
    delete[] wi;
}
