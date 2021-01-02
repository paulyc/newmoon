_: all

all: build

submodules:
	git submodule update --init
.PHONY: submodules

googletest/lib/libgtest.a: submodules
	cd googletest && cmake . && make -j

googletest: googletest/lib/libgtest.a

cmake: CMakeLists.txt
	mkdir -p build
	cd build && cmake ..

build: cmake
	cd build && make -j
.PHONY: build

build/test/test: build test/*
build/src/newmoon: build src/*

run: ephem build/src/newmoon
	build/src/newmoon
.PHONY: run

test: build/test/test
	build/test/test
.PHONY: test

clean:
	rm -rf build
.PHONY: clean

buildtest: build test
.PHONY: buildtest

cleanconf:
	rm -rf build
.PHONY: cleanconf

ephem:
	make -C ephem

ephem-clean:
	make -C ephem clean

distclean: clean ephem-clean
.PHONY: ephem ephem-clean distclean
