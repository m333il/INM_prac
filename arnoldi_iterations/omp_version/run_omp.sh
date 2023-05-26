#!/bin/bash

mkdir -p build
cd build 

cmake ..
make

./arnoldi 1000 100 1
./arnoldi 2000 100 1
./arnoldi 4000 100 1
./arnoldi 6000 100 1
./arnoldi 8000 100 1
./arnoldi 10000 100 1
./arnoldi 12000 100 1
./arnoldi 14000 100 1
./arnoldi 16000 100 1

./arnoldi 1000 100 2
./arnoldi 2000 100 2
./arnoldi 4000 100 2
./arnoldi 6000 100 2
./arnoldi 8000 100 2
./arnoldi 10000 100 2
./arnoldi 12000 100 2
./arnoldi 14000 100 2
./arnoldi 16000 100 2

./arnoldi 1000 100 4
./arnoldi 2000 100 4
./arnoldi 4000 100 4
./arnoldi 6000 100 4
./arnoldi 8000 100 4
./arnoldi 10000 100 4
./arnoldi 12000 100 4
./arnoldi 14000 100 4
./arnoldi 16000 100 4

./arnoldi 1000 100 6
./arnoldi 2000 100 6
./arnoldi 4000 100 6
./arnoldi 6000 100 6
./arnoldi 8000 100 6
./arnoldi 10000 100 6
./arnoldi 12000 100 6
./arnoldi 14000 100 6
./arnoldi 16000 100 6

cd ..
rm -rf build
