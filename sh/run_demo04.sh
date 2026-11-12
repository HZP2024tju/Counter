# !/bin/bash
#  默认参数     我了个现代化脚本 
N=${1:-200}
A=${2:-1}
B=${3:-2}
CHOOSE=${4:-4}


cd ..
mkdir build 
cd build
cmake -DCMAKE_CXX_FLAGS="-fsanitize=address -g -O0" .. #好舒服啊
make -j10
#./mytest
./demo04 "$N" "$A" "$B"  "$CHOOSE" 
cd ..