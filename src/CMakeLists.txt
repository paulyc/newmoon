
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
	ephemshelper.hpp
	quadmath.h
)

target_link_libraries(newmoon -lquadmath)
