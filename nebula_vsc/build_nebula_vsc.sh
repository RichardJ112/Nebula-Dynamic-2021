#!/bin/sh

#Shell script to build nebula_vsc

clear
mkdir -p build
cd build
cmake ..
make
