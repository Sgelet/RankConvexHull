cmake_minimum_required(VERSION 3.24)
project(DynamicConvexHull)

set(CMAKE_CXX_STANDARD 17)

find_package(CGAL)

add_executable(DynamicConvexHull main.cpp ConvexHull.h)

target_link_libraries(DynamicConvexHull CGAL::CGAL)