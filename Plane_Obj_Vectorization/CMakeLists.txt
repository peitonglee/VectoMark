cmake_minimum_required(VERSION 3.10)

# 设置项目名称和版本
project(Plane_Obj_Vectorization VERSION 1.0)


set(CMAKE_BUILD_TYPE "Debug")
# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
find_package(matplotlib_cpp REQUIRED)

# 添加可执行文件
add_executable(Plane_Obj_Vectorization main.cpp ldbscan.cpp Plane_Obj_Vectorization.cpp)
target_link_libraries(Plane_Obj_Vectorization PUBLIC matplotlib_cpp::matplotlib_cpp)
