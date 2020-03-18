cmake_minimum_required(VERSION 3.10)

project(newmoon)
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/lib/gcc/x86_64-pc-linux-gnu/9.2.1/include)
link_directories(extern/lib)

add_subdirectory(extern)
add_subdirectory(src)
add_subdirectory(test)
