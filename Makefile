CXX ?= g++

newmoon: get_bin.h jpleph.cpp jpleph.h jpl_int.h NewMoon.cpp  NewMoon.hpp
	$(CXX) -o newmoon jpleph.cpp NewMoon.cpp

all: newmoon ephems

ephems:
	mkdir -p ephems
	wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/linux_p1550p2650.430
	wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430t/linux_p1550p2650.430t
	wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de431/lnxm13000p17000.431

clean:
	rm -f newmoon
.PHONY: clean

distclean: clean
	rm -rf ephems
.PHONY: distclean
