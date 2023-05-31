cmake_minimum_required(VERSION 3.10)
project(lcde LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

add_library(pnnls OBJECT "include/lcde/pnnls.f")
file(GLOB INCLUDE_H "include/lcde/*.h")
set(TEST_FILES test.cc)

add_executable(test ${INCLUDE_H} ${TEST_FILES})

target_link_libraries(test PRIVATE pnnls)