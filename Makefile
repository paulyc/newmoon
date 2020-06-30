_: all

all: build

submodules:
	git submodule update --init
.PHONY: submodules

googletest/lib/libgtest.a: submodules
	cd googletest && cmake . && make -j

googletest: googletest/lib/libgtest.a

build: CMakeLists.txt googletest
	rm -rf build
	mkdir -p build
	cd build && cmake .. && make -j

build/test/test: build test/*
build/src/newmoon: build src/*

run: build/src/newmoon
	build/src/newmoon
.PHONY: run

test: build/test/test
	build/test/test
.PHONY: test

clean:
	rm -rf build
	rm -rf googletest
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

ephem-clean:
	rm -rf ephem

distclean: clean ephem-clean
