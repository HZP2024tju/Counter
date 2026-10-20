# !/bin/bash
#  默认参数     我了个现代化脚本 
N=${1:-30}
A=${2:--20.4}
B=${3:-20.4}
C=${4:-4}
D=${5:-2}
E=${6:-2}
F=${7:-1}
CHOOSE=${8:-0}
DEVIDE=${9:-6}
UP=${10:-4} #多项式次数
RAND_SET=${11:-0.0}

cd ..
mkdir build 
cd build
cmake -DCMAKE_CXX_FLAGS="-fsanitize=address -g -O0" .. #好舒服啊
make -j10
#./mytest
./demo03 "$N" "$A" "$B" "$C" "$D" "$E" "$F" "$CHOOSE" "$DEVIDE" "$UP" "$RAND_SET"
cd ..
