cmake_minimum_required(VERSION 3.10)
project(lcde LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

add_library(pnnls OBJECT "include/lcde/pnnls.f")

add_executable(build build.cc)
target_link_libraries(build PRIVATE pnnls)

add_executable(error error.cc)
target_link_libraries(error PRIVATE pnnls)