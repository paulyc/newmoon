cmake_minimum_required(VERSION 3.10)

project(newmoon)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(USE_FLOAT128 "-DUSE_FLOAT128=0")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${USE_FLOAT128}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USE_FLOAT128}")

unset(INSTALL_GTEST)
add_subdirectory(googletest)
include_directories(googletest/include)
link_directories(googletest/lib)
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(test)
