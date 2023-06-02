cmake_minimum_required(VERSION 3.10)
project(lcde LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

add_library(pnnls OBJECT "include/lcde/pnnls.f")
file(GLOB INCLUDE_H "include/lcde/*.h")
set(OBSOLETE_TEST_FILES obsolete.cc)
set(DEBUG_FILES debug.cc)

add_executable(obsolete ${INCLUDE_H} ${OBSOLETE_TEST_FILES})
add_executable(debug ${INCLUDE_H} ${DEBUG_FILES})

target_link_libraries(obsolete PRIVATE pnnls)
target_link_libraries(debug PRIVATE pnnls)