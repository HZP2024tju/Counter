# !/bin/bash
cd ..
mkdir build 
cd build
cmake ..
make -j7
#./mytest
./demo01 30  -10.1 10.1 4 2 2 11 0
cd ..
