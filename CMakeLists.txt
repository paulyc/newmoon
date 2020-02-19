cmake_minimum_required(VERSION 3.10)

project(newmoon)

add_executable(newmoon NewMoon.cpp NewMoon.hpp get_bin.h jpleph.cpp jpleph.h jpl_int.h)
