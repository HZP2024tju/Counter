cmake_minimum_required(VERSION 3.10)
project(counter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

add_subdirectory(thirdparty/matplotlib-cpp) #这只是一个绘图库
include_directories(
  include
)

#我想试试看使用openmp加速
add_executable(demo01 src/demo01.cpp)
target_link_libraries(demo01 PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen matplotlib_cpp Python3::Python)
target_include_directories(demo01 PRIVATE
    ${Python3_INCLUDE_DIRS}
)


#我想试试看使用openmp加速
add_executable(demo02 src/demo02.cpp)
target_link_libraries(demo02 PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen matplotlib_cpp Python3::Python)
target_include_directories(demo02 PRIVATE
    ${Python3_INCLUDE_DIRS}
)

add_executable(demo03 src/demo03.cpp)
target_link_libraries(demo03 PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen matplotlib_cpp Python3::Python)
target_include_directories(demo03 PRIVATE
    ${Python3_INCLUDE_DIRS}
)

add_executable(mytest src/mytest.cpp)
target_link_libraries(mytest PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen matplotlib_cpp Python3::Python)
target_include_directories(mytest PRIVATE
    ${Python3_INCLUDE_DIRS}
)