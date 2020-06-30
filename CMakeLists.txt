cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER /usr/bin/clang)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)

project(newmoon)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

set(USE_FLOAT128 1)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_FLOAT128=${USE_FLOAT128}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FLOAT128=${USE_FLOAT128}")

unset(INSTALL_GTEST)
add_subdirectory(googletest)

include_directories(googletest/include)
link_directories(googletest/lib)

add_subdirectory(src)
add_subdirectory(test)
