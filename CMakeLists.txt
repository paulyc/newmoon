cmake_minimum_required(VERSION 3.10)

project(newmoon)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

#unset(INSTALL_GTEST)
#add_subdirectory(googletest)

include_directories(include)
link_directories(lib)

add_subdirectory(src)
add_subdirectory(test)
