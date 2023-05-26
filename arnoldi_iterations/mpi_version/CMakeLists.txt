cmake_minimum_required(VERSION 3.10)
project(arnoldi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto -march=native -llapacke")

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)

add_executable(arnoldi main.cpp)

target_link_libraries(arnoldi MPI::MPI_CXX)