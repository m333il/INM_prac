cmake_minimum_required(VERSION 3.10)
project(arnoldi)

set(GCC_COVERAGE_COMPILE_FLAGS "-Ofast -flto -march=native -llapacke")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(BLAS REQUIRED)

add_executable(arnoldi main.cpp)