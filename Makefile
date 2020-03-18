all: build

configure: CMakeLists.txt
	mkdir -p build
	cd build && cmake ..

build: configure
	cd build && make -j
.PHONY: build

run: build
	build/src/newmoon
.PHONY: run

test: build
	build/test/test
.PHONY: test

clean:
	cd build && make clean
.PHONY: clean

buildtest: build test
.PHONY: buildtest

cleanconf:
	rm -rf build
.PHONY: cleanconf

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
