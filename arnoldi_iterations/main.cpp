#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

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

double* matvec(double *A, double *x, double *res, int n) {
    for (unsigned i = 0; i < n; ++i) {
        res[i] = dot(A + i * n, x, n);
    }
    return res;
}

void vec_sub(double *v1, double *v2, double alpha, int n) {
    for (int i = 0; i < n; ++i) {
        v1[i] -= alpha * v2[i];
    }
}

void arnoldi_iteration(double *A, double *b, int n, int k, double *Q, double *h) {
    //A - n x n array
    //b - initial vector, length n
    //k - amount of iterations

    double eps = 1e-12;
    double *v = new double[n];
    normalize(b, n);
    for (int i = 0; i < n; ++i) {
        Q[i] = b[i];
    }
    
    for (int i = 1; i <= k; ++i) {
        matvec(A, Q + (i - 1) * n, v, n);
        for (int j = 0; j < i; ++j) {
            h[j * k + i - 1] = dot(Q + j * n, v, n);
            vec_sub(v, Q + j * n, h[j * k + i - 1], n);
        }
        h[i * k + i - 1] = norm(v, n);
        if (h[i * k + i - 1] > eps) {
            //Q[:, k] = v / h[k, k - 1]
            for (int j = 0; j < n; ++j) {
                Q[i * n + j] = v[j] / h[i * k + i - 1];
            }
        } else {
            break;
        }
    }
    delete[] v;
}

double drand(double a, double b) {
    double f = (double)rand() / RAND_MAX;
    return a + f * (b - a);
}

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    
    double a = -10.0, b = 10.0;
    double *A = new double[n * n];
    for (int i = 0; i < n * n; ++i) {
        A[i] = drand(a, b);
    }
    double *v1 = new double[n];
    for (int i = 0; i < n; ++i) {
        v1[i] = drand(a, b);
    }
    double *h = new double[(k + 1) * k];
    double *Q = new double[(k + 1) * n];
    for (int i = 0; i < (k + 1) * k; ++i)
        h[i] = 0.0;
    for (int i = 0; i < (k + 1) * n; ++i)
        Q[i] = 0.0;
    arnoldi_iteration(A, v1, n, k, Q, h);
    std::ofstream A_file, h_file;
    A_file.open("A.txt");
    h_file.open("h.txt");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_file << A[i * n + j] << " ";
        }
        A_file << std::endl;
    }
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            h_file << h[i * k + j] << " ";
        }
        h_file << std::endl;
    }
    A_file.close();
    h_file.close();
    delete[] Q;
    delete[] v1;
    delete[] A;
    delete[] h;
}