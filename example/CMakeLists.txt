# ----- Find dependencies -----

# popl
set(popl_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/3rd/popl/include)

if(NOT USE_PANGOLIN_VIEWER)
    message(FATAL_ERROR "run_slam_plane_vector requires USE_PANGOLIN_VIEWER")
endif()

# ----- Build the vectoMark executable -----
add_executable(run_slam_plane_vector run_slam_plane_vector.cc util/planeSeg_util.cc)

set_target_properties(run_slam_plane_vector PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${PROJECT_BINARY_DIR}")

target_compile_definitions(run_slam_plane_vector PRIVATE USE_PANGOLIN_VIEWER)
target_link_libraries(run_slam_plane_vector PRIVATE pangolin_viewer)

target_link_libraries(run_slam_plane_vector
    PRIVATE
    ${PROJECT_NAME}
    opencv_imgcodecs
    opencv_videoio)

target_include_directories(run_slam_plane_vector
    PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/3rd/popl/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/3rd/spdlog/include>)
