
CMAKE := cmake
CMAKE_CMD := $(CMAKE) -DCMAKE_C_COMPILER=/usr/local/bin/gcc-10 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-10
ROOT := $(PWD)

_: all

all: build

submodules: googletest
	git submodule update --init

lib/libgtest.a: submodules
	cd googletest && rm -rf build && mkdir build && cd build && $(CMAKE_CMD) -DCMAKE_INSTALL_PREFIX=$(ROOT) .. && make -j install

googletest: lib/libgtest.a

cmake: googletest CMakeLists.txt
	mkdir -p build
	cd build && $(CMAKE_CMD) ..
.PHONY: cmake

build: cmake
	cd build && make -j
.PHONY: build

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
	make -C googletest clean
.PHONY: clean

buildtest: build test
.PHONY: buildtest

cleanconf:
	rm -rf build
.PHONY: cleanconf

ephem: ephem/lnxm13000p17000.431
	mkdir -p ephem
	cd ephem && wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de431/lnxm13000p17000.431

ephem-clean:
	rm -rf ephem

distclean: clean ephem-clean
