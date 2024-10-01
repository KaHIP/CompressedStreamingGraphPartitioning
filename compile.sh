#!/bin/bash

NCORES=12
unamestr=`uname`

rm -rf deploy
rm -rf build
mkdir build
cd build 
cmake -DCMAKE_CXX_COMPILER=g++ ../
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg ../
make -j $NCORES
cd ..

mkdir deploy
cp ./build/stream_cpi deploy/
cp ./build/stream_cpi_generated deploy/

rm -rf build
