# !/bin/bash
#  默认参数     我了个现代化脚本 
cd ..
mkdir build 
cd build
cmake -DCMAKE_CXX_FLAGS="-fsanitize=address -g -O0" .. #好舒服啊
make -j10
#       max_step  err w
./demo06 30000 1e-8
cd ..