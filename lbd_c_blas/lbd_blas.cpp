#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

#include <cblas.h>

unsigned n;

double* transpose(double **A) {
    double *AT = new double[n * n];
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            AT[i * n + j] = A[j][i];
    }
    return AT;
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

void reorthog(double *v, double **M, int it) {
    for (int i = 1; i < it; ++i) {
        double a = cblas_ddot(n, M[i], 1, v, 1);
	//v = -a * M[i] + v	
	cblas_daxpy(n, -a, M[i], 1, v, 1);
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

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (unsigned k = 1; k <= n; ++k) {
        //r = A * V[k] - betta[k - 1] * U[k - 1]
        double *r = new double[n];
        cblas_dcopy(n, U[k - 1], 1, r, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, V[k], 1, -betta[k - 1], r, 1);
        reorthog(r, U, k);
        U[k] = r;
        alpha[k] = cblas_dnrm2(n, U[k], 1);
	//U[k] /= norm(U[k])
        cblas_dscal(n, 1.0 / alpha[k], U[k], 1);
        
	//p = A^T * U[k] - alpha[k] * V[k]
        double *p = new double[n];
        cblas_dcopy(n, V[k], 1, p, 1);
	cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0, A, n, U[k], 1, -alpha[k], p, 1);
        reorthog(p, V, k + 1);
        V[k + 1] = p;
        betta[k] = cblas_dnrm2(n, V[k + 1], 1);
        //V[k + 1] /= norm(V[k + 1])
	cblas_dscal(n, 1.0 / betta[k], V[k + 1], 1);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    return Result(alpha, betta, U, V, A, std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

void check(double *A, double *U, double *V, double *B) {
    //orthogonality
    double *UUT = new double[n * n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, U, n, U, n, 0.0, UUT, n);
    for (int i = 0; i < n; ++i)
        UUT[i * n + i] -= 1.0;
    double sum2 = cblas_dnrm2(n * n, UUT, 1);
    std::cout << "norm(I - UU^T) = " << sum2 << std::endl;

    double *VVT = new double[n * n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, V, n, V, n, 0.0, VVT, n);
    for (int i = 0; i < n; ++i)
        VVT[i * n + i] -= 1.0;
    sum2 = cblas_dnrm2(n * n, VVT, 1);
    std::cout << "norm(I - VV^T) = " << sum2 << std::endl;

    //norm(AV - UB)
    double *AV = new double[n * n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, V, n, 0.0, AV, n);
    double *UB = new double[n * n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, U, n, B, n, 0.0, UB, n);
    cblas_daxpy(n * n, -1.0, AV, 1, UB, 1);
    sum2 = cblas_dnrm2(n * n, UB, 1);
    std::cout << "norm(AV - UB) = " << sum2 << std::endl;
    delete[] AV;
    delete[] UB;
    delete[] UUT;
    delete[] VVT;
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
    double pnorm = cblas_dnrm2(n, p, 1);
    cblas_dscal(n, 1.0 / pnorm, p, 1);
    auto Res = lbd(A, p);
    std::cout << "Time = " << Res.time / 1000000000.0 << " sec" << std::endl;
    check(Res.A, Res.U, Res.V, Res.B);
    
    std::ofstream foutTime;
    foutTime.open("Time.txt");
    foutTime << Res.time;

    foutTime.close();
}
