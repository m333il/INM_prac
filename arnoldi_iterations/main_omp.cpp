#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <lapacke.h>

int num_threads;

double drand(double a, double b) {
    double f = (double)rand() / RAND_MAX;
    return a + f * (b - a);
}

class Matrix
{
private:
    const uint32_t n;
    double *A;
    const double a, b;
public:
    Matrix(uint32_t size, double a = -10.0, double b = 10.0) : 
        n(size), a(a), b(b)
    {
        A = new double[n * n];
        for (uint32_t i = 0; i < n * n; ++i) {
            A[i] = drand(a, b);
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

    void matvec(double *x, double *res) const {
        uint32_t i, j;
        #pragma omp parallel for num_threads(num_threads)
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

double norm(double *vec, int n) {
    double res = 0.0;
    for (unsigned i = 0; i < n; ++i)
        res += vec[i] * vec[i];
    return sqrt(res);
}

void normalize(double *vec, int n) {
    double vnorm = norm(vec, n);
    for (int i = 0; i < n; ++i) {
        vec[i] /= vnorm;
    }
}

double dot(double *v1, double *v2, int n) {
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

void vec_sub(double *v1, double *v2, double alpha, int n) {
    
    for (int i = 0; i < n; ++i) {
        v1[i] -= alpha * v2[i];
    }
}

void arnoldi_iteration(Matrix &A, int k, double *Q, double *h) {
    //A - n x n array
    //b - initial vector, length n
    //k - amount of iterations

    const uint32_t n = A.order();
    double eps = 1e-12;
    double *v = new double[n];

    for (int i = 0; i < n; ++i) {
        Q[i] = drand(-10.0, 10.0);
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
    num_threads = atoi(argv[3]);
    Matrix A(n);

    double *Q = new double[(k + 1) * n];
    double *h = new double[(k + 1) * k];
    std::fill(Q, Q + (k + 1) * n, 0);
    std::fill(h, h + (k + 1) * k, 0);
    double *wr = new double[k];
    double *wi = new double[k];
    double *Z = new double[k * k];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    arnoldi_iteration(A, k, Q, h);

    //std::ofstream h_file;
    //A.save("A.txt");
    //h_file.open("h.txt");
    /*
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            h_file << h[i * k + j] << " ";
        }
        h_file << std::endl;
    }
    h_file.close();
    */
    
    int info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'E', 'N', k, 1, k, h, k, wr, wi, Z, k);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference for n = " << n << " and num_threads = " << num_threads << " is  " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
    /*
    std::ofstream wr_txt, wi_txt;
    wr_txt.open("wr.txt");
    wi_txt.open("wi.txt");
    for (int i = 0; i < k; ++i) {
        wr_txt << wr[i] << " ";
        wi_txt << wi[i] << " ";
    }
    wr_txt.close();
    wi_txt.close();
    */
    delete[] Z;
    delete[] Q;
    delete[] h;
    delete[] wr;
    delete[] wi;
}
