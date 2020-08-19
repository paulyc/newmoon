cmake_minimum_required(VERSION 3.10)

project(newmoon)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_C_COMPILER /usr/bin/clang)
#set(CMAKE_CXX_COMPILER /usr/bin/clang++)

#unset(INSTALL_GTEST)
#add_subdirectory(googletest)

include_directories(googletest/include)
#link_directories(googletest/lib)

add_subdirectory(src)
add_subdirectory(test)
