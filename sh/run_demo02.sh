# !/bin/bash
cd ..
mkdir build 
cd build
cmake ..
make -j10
#./mytest
./demo02 20  -1.1 1.1 4 2 2 11 2
cd ..
