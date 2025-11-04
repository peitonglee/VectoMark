cmake_minimum_required(VERSION 3.10)

# 设置项目名称和版本
project(Path_Planning_Base_Lines VERSION 1.0)


set(CMAKE_BUILD_TYPE "Debug")
# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
find_package(matplotlib_cpp REQUIRED)

# 添加可执行文件
add_executable(Path_Planning_Base_Lines main.cpp)
target_link_libraries(Path_Planning_Base_Lines PUBLIC matplotlib_cpp::matplotlib_cpp)