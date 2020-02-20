cmake_minimum_required(VERSION 3.10)

project(newmoon)

set(CMAKE_CXX_STANDARD 17)

add_executable(newmoon
	newmoon.cpp
	newmoon.hpp
	get_bin.h
	jpleph.cpp
	jpleph.h
	jpl_int.h
	jd_clock.hpp
	jd_clock.cpp
	lalgebra.hpp 
	ephemshelper.hpp
)
