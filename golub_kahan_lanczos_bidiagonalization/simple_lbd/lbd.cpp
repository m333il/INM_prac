#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

unsigned n;

double norm(double *vec) {
    double res = 0.0;
    for (unsigned i = 0; i < n; ++i)
        res += vec[i] * vec[i];
    return sqrt(res);
}

double* transpose(double *A) {
    double *AT = new double[n * n];
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            AT[j * n + i] = A[i * n + j];
    return AT;
}

double* transpose(double **A) {
    double *AT = new double[n * n];
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            AT[i * n + j] = A[j][i];
    }
    return AT;
}

double dot(double *v1, double *v2) {
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

struct Result
{
    double *A;
    double *U;
    double *V;
    double *B;
    long long unsigned time;

    Result(double *alpha, double *betta, double **_U, double **_V, double *_A, long long unsigned _time) {
        double *new_alpha = new double[n];
        double *new_betta = new double[n];
        for (unsigned i = 0; i < n; ++i) {
            new_alpha[i] = alpha[i + 1];
            new_betta[i] = betta[i + 1];
        }
            
        delete[] alpha;
        delete[] betta;

        double **new_U = new double*[n];
        double **new_V = new double*[n];
        for (unsigned i = 0; i < n; ++i) {
            new_U[i] = _U[i + 1];
            new_V[i] = _V[i + 1];
        }

        A = _A;
        U = transpose(new_U);
        V = transpose(new_V);

        delete[] new_U;
        delete[] new_V;
        for (int i = 0; i < n + 1; ++i)
            delete[] _U[i];
        delete[] _U;
        for (int i = 0; i < n + 2; ++i)
            delete[] _V[i];
        delete[] _V;
        time = _time;

        B = new double[n * n];
        for (unsigned i = 0; i < n * n; ++i)
            B[i] = 0.0;
        for (unsigned i = 0; i < n - 1; ++i) {
            B[i * n + i] = new_alpha[i];
            B[i * n + i + 1] = new_betta[i];
        }
        B[n * n - 1] = new_alpha[n - 1];

        delete[] new_alpha;
        delete[] new_betta;
    }

    ~Result() {
        delete[] A;
        delete[] B;
        delete[] U;
        delete[] V;
    }
};

double* multiply_substract(double *A, double *Vk, double *Uk, double betta) {
    double *res = new double[n]();
    for (unsigned i = 0; i < n; ++i) {
        res[i] = dot(A + i * n, Vk);
        //for (unsigned j = 0; j < n; ++j) {
          //  res[i] += A[i * n + j] * Vk[j];
        //}
    }
    for (unsigned i = 0; i < n; ++i) {
        res[i] -= betta * Uk[i];
    }
    return res;
}

void reorthog(double *v, double **M, int it) {
    for (int i = 1; i < it; ++i) {
        double a = dot(M[i], v);
        for (int j = 0; j < n; ++j) {
            v[j] -= a * M[i][j];
        }
    }
}

Result lbd(double* A, double* v1) {
    double *betta = new double[n + 1]();
    double *alpha = new double[n + 1]();

    double **U = new double*[n + 1];
    double **V = new double*[n + 2];

    V[1] = v1;
    U[0] = new double[n]();
    V[0] = new double[n]();

    auto AT = transpose(A);
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (unsigned k = 1; k <= n; ++k) {
        auto r = multiply_substract(A, V[k], U[k - 1], betta[k - 1]);
        reorthog(r, U, k);
        U[k] = r;
        alpha[k] = norm(U[k]);

        for (unsigned i = 0; i < n; ++i)
            U[k][i] /= alpha[k];

        auto p = multiply_substract(AT, U[k], V[k], alpha[k]);
        reorthog(p, V, k + 1);
        V[k + 1] = p;
        betta[k] = norm(V[k + 1]);

        for (unsigned i = 0; i < n; ++i) {
            V[k + 1][i] /= betta[k];
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    delete[] AT;
    return Result(alpha, betta, U, V, A, std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

double* matmul(double *A, double *B) {
    double *C = new double[n * n]();
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    return C;
}

void check(double *A, double *U, double *V, double *B) {
    //orthogonality
    double sum2 = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double curr_el = 0.0;
            for (int k = 0; k < n; ++k) {
                curr_el += U[i * n + k] * U[j * n + k];
            }
            if (i == j) {
                curr_el -= 1.0;
            }
            sum2 += curr_el * curr_el;
        }
    }
    std::cout << "norm(I - UU^T) = " << sqrt(sum2) << std::endl;
    sum2 = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double curr_el = 0.0;
            for (int k = 0; k < n; ++k) {
                curr_el += V[i * n + k] * V[j * n + k];
            }
            if (i == j) {
                curr_el -= 1.0;
            }
            sum2 += curr_el * curr_el;
        }
    }
    std::cout << "norm(I - VV^T) = " << sqrt(sum2) << std::endl;

    //norm(AV - UB)
    double *AV = matmul(A, V);
    double *UB = matmul(U, B);
    sum2 = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double curr = AV[i * n + j] - UB[i * n + j];
            sum2 += curr * curr;
        }
    }
    std::cout << "norm(AV - UB) = " << sqrt(sum2) << std::endl;
    delete[] AV;
    delete[] UB;
}

double drand(double a, double b)
{
    double f = (double)rand() / RAND_MAX;
    return a + f * (b - a);
}

int main(int argc, char *argv[]) {
    n = atoi(argv[1]);
    double a = -10.0, b = 10.0;
    double *A = new double[n * n];
    for (int i = 0; i < n * n; ++i) {
        A[i] = drand(a, b);
    }

    double *p = new double[n];
    for (int i = 0; i < n; ++i) {
        p[i] = drand(a, b);
    }
    double pnorm = norm(p);
    for (int i = 0; i < n; ++i) {
        p[i] /= pnorm;
    }
    auto Res = lbd(A, p);
    std::cout << "Time = " << Res.time / 1000000000.0 << " sec" << std::endl;
    check(Res.A, Res.U, Res.V, Res.B);
    
    std::ofstream foutTime;
    foutTime.open("Time.txt");
    foutTime << Res.time;

    foutTime.close();
}
