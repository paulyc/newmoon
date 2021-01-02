
project(newmoon_executable)

add_executable(newmoon
	newmoon.cpp
	get_bin.h
	jpleph.cpp
	jpleph.h
	jpl_int.h
	jd_clock.hpp
	lalgebra.hpp
	calculus.hpp
	calculus.cpp
	ephemshelper.hpp
	quadmath.h
	tetrabiblos.hpp
	tetrabiblos.cpp
	astro.hpp
	astro.cpp
)

target_link_libraries(newmoon -lquadmath)
