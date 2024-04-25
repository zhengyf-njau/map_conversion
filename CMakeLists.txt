cmake_minimum_required(VERSION 3.0)
project(ModifyPointCloud)

find_package(PCL 1.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(modify_point_cloud modify_point_cloud.cpp)
target_link_libraries(modify_point_cloud ${PCL_LIBRARIES})