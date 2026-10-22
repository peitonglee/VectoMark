cmake_minimum_required(VERSION 3.8)
project(PLPSLAM LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 OLD)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# ----- Set build type -----
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ----- Set options for debugging -----
set(USE_CCACHE ON CACHE BOOL "Use ccache to accelerate build")
find_program(CCACHE_EXE ccache)

if(USE_CCACHE AND CCACHE_EXE)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_EXE}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_EXE}")
    message(STATUS "ccache: Enabled")
else()
    set(USE_CCACHE OFF)
    message(STATUS "ccache: Disabled")
endif()

# ----- Set compiler options -----
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Og")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(BUILD_WITH_MARCH_NATIVE OFF CACHE BOOL "Enable architecture-aware optimization (-march=native)")

if(BUILD_WITH_MARCH_NATIVE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -march=native")
    message(STATUS "Architecture-aware optimization (-march=native): ENABLED")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native")
    message(STATUS "Architecture-aware optimization (-march=native): DISABLED")
endif()

# ----- Eigen alignment workaround for g2o -----
# Must match the macro used to build the system g2o libraries (/usr/local/lib/libg2o_*.so).
add_definitions(-DEIGEN_DONT_ALIGN_STATICALLY)

# ----- Find dependencies -----

# Threads
find_package(Threads REQUIRED)

# Eigen
find_package(Eigen3 3.3 REQUIRED)

# yaml-cpp
find_package(yaml-cpp REQUIRED)

# OpenCV
find_package(OpenCV 3.3.1 QUIET
    COMPONENTS
    core imgcodecs videoio features2d calib3d highgui)

if(NOT OpenCV_FOUND)
    find_package(OpenCV 4.0 REQUIRED
        COMPONENTS
        core imgcodecs videoio features2d calib3d highgui)

    if(NOT OpenCV_FOUND)
        message(FATAL_ERROR "OpenCV >= 3.3.1 not found")
    endif()
endif()

message(STATUS "Use OpenCV ${OpenCV_VERSION}")

# ----- Install configuration -----
# (kept so that the install() rules inside src/ and its subdirectories stay valid)
include(CMakePackageConfigHelpers)

set(PLPSLAM_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(PLPSLAM_VERSION_CONFIG ${PLPSLAM_GENERATED_DIR}/${PROJECT_NAME}ConfigVersion.cmake)
set(PLPSLAM_PROJECT_CONFIG ${PLPSLAM_GENERATED_DIR}/${PROJECT_NAME}Config.cmake)
set(PLPSLAM_TARGETS_EXPORT_NAME ${PROJECT_NAME}Targets)
set(PLPSLAM_CONFIG_INSTALL_DIR lib/cmake/${PROJECT_NAME})
set(PLPSLAM_NAMESPACE "${PROJECT_NAME}::")
set(PLPSLAM_VERSION 0.0.1)

write_basic_package_version_file(${PLPSLAM_VERSION_CONFIG}
    VERSION ${PLPSLAM_VERSION}
    COMPATIBILITY SameMajorVersion)

configure_file(${PROJECT_SOURCE_DIR}/cmake/Config.cmake.in ${PLPSLAM_PROJECT_CONFIG} @ONLY)

install(FILES ${PLPSLAM_VERSION_CONFIG} ${PLPSLAM_PROJECT_CONFIG}
    DESTINATION ${PLPSLAM_CONFIG_INSTALL_DIR})
install(EXPORT ${PLPSLAM_TARGETS_EXPORT_NAME}
    NAMESPACE ${PLPSLAM_NAMESPACE}
    DESTINATION ${PLPSLAM_CONFIG_INSTALL_DIR})

set(RUNTIME_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
set(LIBRARY_DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
set(ARCHIVE_DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
set(INCLUDES_DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

# ----- Build selection -----
set(USE_PANGOLIN_VIEWER ON CACHE BOOL "Enable Pangolin Viewer")

add_subdirectory(src)
add_subdirectory(example)
