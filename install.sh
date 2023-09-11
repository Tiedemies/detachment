#!/bin/bash
cd src 
cmake -DCMAKE_BUILD_TYPE=Release
make clean
make install
cd ..

