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
            AT[i + n * j] = A[j + i * n];
    return AT;
}

double** transpose(double **A) {
    double **AT = new double*[n];
    for (size_t i = 0; i < n; ++i) {
        AT[i] = new double[n];
        for (size_t j = 0; j < n; ++j)
            AT[i][j] = A[j][i];
    }
    return AT;
}

struct Result
{
    double *A;
    double **U;
    double **V;
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
        delete[] _U;
        delete[] _V;
        time = _time;

        B = new double[n * n];
        for (unsigned i = 0; i < n * n; ++i)
            B[i] = 0.0;
        for (unsigned i = 0; i < n - 1; ++i) {
            B[i * n + i] = new_alpha[i];
            B[i + (i + 1) * n] = new_betta[i];
        }
        B[n - 1 + (n - 1) * n] = new_alpha[n - 1];

        delete[] new_alpha;
        delete[] new_betta;
    }

    ~Result() {
        delete[] A;
        delete[] U;
        delete[] V;
        delete[] B;
    }
};

double* multiply_substract(double *A, double *Vk, double *Uk, double betta) {
    double *res = new double[n];
    for (unsigned j = 0; j < n; ++j) {
        for (unsigned i = 0; i < n; ++i) {
            res[i] += A[i + j * n] * Vk[j];
        }
    }
    for (unsigned i = 0; i < n; ++i) {
        res[i] -= betta * Uk[i];
    }
    return res;
}

Result lbd(double* A, double* v1) {
    double *betta = new double[n + 1];
    double *alpha = new double[n + 1];
    for (unsigned i = 0; i < n + 1; ++i) {
        betta[i] = 0.0;
        alpha[i] = 0.0;
    }

    double **U = new double*[n + 1];
    double **V = new double*[n + 2];

    for (unsigned i = 0; i < n + 1; ++i) {
        U[i] = new double[n];
        V[i] = new double[n];
        for (unsigned j = 0; j < n; ++j) {
            U[i][j] = 0.0;
            V[i][j] = 0.0;
        }
    }
    V[n + 1] = new double[n];
    for (unsigned i = 0; i < n; ++i)
        V[1][i] = v1[i];

    auto AT = transpose(A);
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (unsigned k = 1; k <= n; ++k) {
        U[k] = multiply_substract(A, V[k], U[k - 1], betta[k - 1]);
        alpha[k] = norm(U[k]);

        for (unsigned i = 0; i < n; ++i)
            U[k][i] /= alpha[k];

        V[k + 1] = multiply_substract(AT, U[k], V[k], alpha[k]);
        betta[k] = norm(V[k + 1]);

        for (unsigned i = 0; i < n; ++i) {
            V[k + 1][i] /= betta[k];
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    delete[] AT;
    delete[] v1;
    return Result(alpha, betta, U, V, A, std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

int main(int argc, char *argv[]) {
    n = atoi(argv[1]);

    double *A = new double[n * n];
    std::ifstream fin;
    fin.open("test.txt");
    for (size_t i = 0; i < n * n; ++i) {
        fin >> A[i];
    }
    fin.close();

    double *p = new double[n];
    p[0] = 1.0;
    for (unsigned i = 1; i < n; ++i)
        p[i] = 0.0;
    auto Res = lbd(A, p);

    std::ofstream foutA;
    foutA.open("A.txt");

    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            foutA << Res.A[i + j * n] << " ";
        }
        foutA << std::endl;
    }

    std::ofstream foutU;
    foutU.open("U.txt");

    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            foutU << Res.U[i][j] << " ";
        }
        foutU << std::endl;
    }

    std::ofstream foutV;
    foutV.open("V.txt");
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            foutV << Res.V[i][j] << " ";
        }
        foutV << std::endl;
    }

    std::ofstream foutB;
    foutB.open("B.txt");
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            foutB << Res.B[i + j * n] << " ";
        }
        foutB << std::endl;
    }

    std::ofstream foutTime;
    foutTime.open("Time.txt");
    foutTime << Res.time;

    foutA.close();
    foutB.close();
    foutU.close();
    foutV.close();
    foutTime.close();
}
