cmake_minimum_required(VERSION 3.10)
project(lcde LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

# include_directories(
#   ${CMAKE_SOURCE_DIR}
# )

add_library(pnnls OBJECT "src/pnnls.f")
file(GLOB INCLUDE_H "src/*.h")
set(TEST_FILES test.cc)

add_executable(test ${INCLUDE_H} ${TEST_FILES})

target_link_libraries(test PRIVATE pnnls)