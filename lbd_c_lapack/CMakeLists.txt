project(bidiagonalization)

cmake_minimum_required(VERSION 3.24)
set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapacke")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(BLAS REQUIRED)

add_executable(lbd_lapack lbd_lapack.cpp)
