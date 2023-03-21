#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

#include <lapacke.h>

unsigned n;

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
    double *D = new double[n];
    double *E = new double[n - 1];
    double *tauq = new double[n];
    double *taup = new double[n];
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    int info = LAPACKE_dgebrd(LAPACK_ROW_MAJOR, n, n, A, n, D, E, tauq, taup);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    long long unsigned time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    std::cout << "Time = " << time / 1000000000.0 << " sec" << std::endl;
    
    std::ofstream foutTime;
    foutTime.open("Time.txt");
    foutTime << time;

    foutTime.close();
    delete[] D;
    delete[] E;
    delete[] tauq;
    delete[] taup;
    delete[] A;
}
