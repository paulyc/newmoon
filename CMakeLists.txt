cmake_minimum_required(VERSION 3.10)

project(newmoon)
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(extern)
add_subdirectory(src)
add_subdirectory(test)
