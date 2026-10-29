cmake_minimum_required(VERSION 3.10)
project(counter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # 给 VS Code 用的神器，保留！

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

add_subdirectory(thirdparty/matplotlib-cpp)


# 公共接口库
add_library(common_deps INTERFACE)

target_link_libraries(common_deps INTERFACE
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
    matplotlib_cpp
    Python3::Python
)

target_include_directories(common_deps INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include                 
    ${Python3_INCLUDE_DIRS}                              
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/matplotlib-cpp 
)

# 批量处理普通 demo 
set(NORMAL_TARGETS demo01 demo02 demo03 demo05 mytest)
foreach(tgt ${NORMAL_TARGETS})
    add_executable(${tgt} src/${tgt}.cpp)
    target_link_libraries(${tgt} PRIVATE common_deps)
endforeach()

# 单独处理带内存检查 的 demo04

add_executable(demo04 src/demo04.cpp)
target_link_options(demo04 PRIVATE -fsanitize=address -static-libasan)
target_link_libraries(demo04 PRIVATE common_deps)