cmake_minimum_required(VERSION 3.24)
project(DynamicConvexHull)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_BUILD_TYPE Release)

#find_package(CGAL)

add_executable(DynamicConvexHull main.cpp ConvexHull.h)

#target_link_libraries(DynamicConvexHull CGAL::CGAL)