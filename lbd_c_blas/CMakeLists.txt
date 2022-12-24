project(bidiagonalization)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -lcblas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(BLAS REQUIRED)

add_executable(lbd_blas lbd_blas.cpp)
