#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

double norm(const std::vector<double> &vec) {
    double res = 0.0;
    for (const auto &el : vec) {
        res += el * el;
    }
    return sqrt(res);
}

std::vector<double> operator/(std::vector<double> vec, double r) {
    for (auto &el : vec)
        el /= r;
    return vec;
}

std::vector<double> operator*(std::vector<double> vec, double r) {
    for (auto &el : vec)
        el *= r;
    return vec;
}


std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &A) {
    size_t rows = A.size();
    size_t columns = A[0].size();
    std::vector<std::vector<double>> AT(columns);
    for (size_t i = 0; i < columns; ++i)
        for (size_t j = 0; j < rows; ++j)
            AT[i].push_back(A[j][i]);
    return AT;
}

double operator*(const std::vector<double> &v1, const std::vector<double> &v2) {
    double res = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
        res += v1[i] * v2[i];
    return res;
}

std::vector<double> operator*(const std::vector<std::vector<double>> &A, const std::vector<double> &v) {
    std::vector<double> res(v.size(), 0.0);
    for (size_t i = 0; i < v.size(); ++i)
        res[i] = A[i] * v;
    return res;
}

template<class T, class P>
std::vector<std::vector<double>> operator*(const std::vector<std::vector<T>> &A, const std::vector<std::vector<P>> &B) {
    std::vector<std::vector<T>> res(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        res[i] = std::vector<double>(A.size(), 0.0);
        for (size_t j = 0; j < A.size(); ++j) {
            for (size_t s = 0; s < A.size(); ++s) {
                res[i][j] += A[i][s] * B[s][j];
            }
        }
    }
    return res;
}


std::vector<double> operator-(std::vector<double> v1, const std::vector<double> &v2) {
    for (size_t i = 0; i < v2.size(); ++i)
        v1[i] -= v2[i];
    return v1;
}

template <class T>
void print_vec(const std::vector<T> &vec) {
    for (const auto &el : vec)
        std::cout << el << " ";
    std::cout << std::endl;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &M) {
    for (const auto &el : M)
        os << el << " ";
    return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>> &M) {
    for (const auto &el : M)
        os << el << std::endl;
    return os;
}

struct Result
{
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> U;
    std::vector<std::vector<double>> V;
    std::vector<std::vector<double>> B;
    long long unsigned time;

    Result(const std::vector<double> &alpha, const std::vector<double> &betta, const std::vector<std::vector<double>> &_U, const std::vector<std::vector<double>> &_V, const std::vector<std::vector<double>> &_A, long long unsigned _time) {
        A = _A;
        U = transpose(_U);
        V = transpose(_V);
        time = _time;
        B.resize(A.size());
        for (unsigned i = 0; i < A.size() - 1; ++i) {
            B[i] = std::vector<double>(A.size(), 0.0);
            B[i][i] = alpha[i];
            B[i][i + 1] = betta[i];
        }
        B[A.size() - 1] = std::vector<double>(A.size(), 0.0);
        B[A.size() - 1][A.size() - 1] = alpha[A.size() - 1];
    }

    void check_equation() const {
        std::cout << "Mutiplication result: " << std::endl;
        std::cout << A * V << std::endl;
        std::cout << U * B << std::endl;
    }
};

Result lbd(std::vector<std::vector<double>> A, std::vector<double> v1) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    std::vector<double> betta(A.size() + 1, 0.0);
    betta[0] = 0.0;
    std::vector<double> alpha(A.size() + 1, 0.0);
    std::vector<std::vector<double>> U(A.size() + 1);
    std::vector<std::vector<double>> V(A.size() + 2);
    V[1] = v1;
    
    for (unsigned k = 1; k <= A.size(); ++k) {
        U[k] = A * V[k] - U[k - 1] * betta[k - 1];
        alpha[k] = norm(U[k]);
        U[k] = U[k] / alpha[k];
        V[k + 1] = transpose(A) * U[k] - V[k] * alpha[k];
        betta[k] = norm(V[k + 1]);
        V[k + 1] = V[k + 1] / betta[k];
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    alpha.erase(alpha.begin());
    betta.erase(betta.begin());
    U.erase(U.begin());
    V.erase(V.end() - 1);
    V.erase(V.begin());
    return Result(alpha, betta, U, V, A, std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
}

int main(int argc, char *argv[]) {
    int n = atoi(argv[1]);
    std::ifstream fin;
    fin.open("test.txt");
    std::vector<std::vector<double>> A(n);
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A.size(); ++j) {
            double a;
            fin >> a;
            A[i].push_back(a);
        }
    }
    fin.close();
    std::vector<double> p(n, 0.0);
    p[0] = 1.0;
    auto Res = lbd(A, p);

    std::ofstream fout;
    fout.open("res_test.txt");
    fout << Res.A << std::endl << Res.U << std::endl << Res.V << std::endl << Res.B << std::endl;
    fout << Res.time;
    fout.close();
}
