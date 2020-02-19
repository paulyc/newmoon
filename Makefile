CXX ?= g++

newmoon: get_bin.h jpleph.cpp jpleph.h jpl_int.h NewMoon.cpp NewMoon.hpp
	$(CXX) -o newmoon jpleph.cpp NewMoon.cpp

all: newmoon ephem

ephem:
	mkdir -p ephem
	wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/linux_p1550p2650.430
	wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430t/linux_p1550p2650.430t

ephem431:
	mkdir -p ephem
	wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de431/lnxm13000p17000.431

clean:
	rm -f newmoon
.PHONY: clean

distclean: clean
	rm -rf ephem
.PHONY: distclean
