all: build

build:
	mkdir -p build
	cd build && cmake ..
	cd build && make

run: build
	build/src/newmoon
.PHONY: run

test: build
	build/test/test
.PHONY: test

clean:
	rm -rf build
.PHONY: clean

ephem:
	mkdir -p ephem
	cd ephem && wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/linux_p1550p2650.430
	cd ephem && wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430t/linux_p1550p2650.430t

ephem431:
	mkdir -p ephem
	cd ephem && wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de431/lnxm13000p17000.431

distclean: clean
	rm -rf ephem
.PHONY: distclean
