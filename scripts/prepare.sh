#! /usr/bin/env bash

echo "Compiling LCDE..."

rm -rf build
mkdir -p build
cd build
cmake ..
make -j 8