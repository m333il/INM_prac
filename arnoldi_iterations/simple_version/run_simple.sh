#!/bin/bash

mkdir -p build
cd build 

cmake ..
make

./arnoldi 1000 100
./arnoldi 2000 100
./arnoldi 4000 100
./arnoldi 6000 100
./arnoldi 8000 100
./arnoldi 10000 100
./arnoldi 12000 100
./arnoldi 14000 100
./arnoldi 16000 100

cd ..
rm -rf build
