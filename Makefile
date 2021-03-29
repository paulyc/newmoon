_: build

all: ephem build

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

ephem:
	make -C ephem

ephem-clean:
	make -C ephem clean

distclean: clean ephem-clean
.PHONY: ephem ephem-clean distclean
