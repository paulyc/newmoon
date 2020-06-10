

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

target_include_directories(newmoon PUBLIC ../include)

target_link_libraries(newmoon -lquadmath)
