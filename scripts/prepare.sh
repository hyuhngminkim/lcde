#! /usr/bin/env bash

echo "Compiling LCDE..."

rm -rf build
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 8