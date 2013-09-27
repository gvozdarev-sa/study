

LIB=-lIL


CFLAGS=-fopenmp
CXXFLAGS=${CFLAGS}
CXX=g++

all: dirs z1

obj/z1.o: src/z1.cpp
	${CXX} -c ${CXXFLAGS} src/z1.cpp -o obj/z1.o

obj/image.o: src/Image.cpp src/Image.hpp
	${CXX} -c ${CXXFLAGS} src/Image.cpp -o obj/image.o

z1: obj/z1.o obj/image.o
	${CXX}  ${CXXFLAGS} obj/z1.o obj/image.o -o bin/z1 ${LIB}

dirs:
	mkdir -p bin
	mkdir -p obj


devil:
	wget http://downloads.sourceforge.net/openil/DevIL-1.7.8.tar.gz
	tar -xvf DevIL-1.7.8.tar.gz
	cd  devil-1.7.8 && ./configure
	cd  devil-1.7.8 && make
	cd  devil-1.7.8 && make install DESTDIR=`pwd`/..
	cp -rf usr/local/include/* include
	cp -rf usr/local/lib/* lib
	rm -rf devil-1.7.8
	rm -rf usr
	rm DevIL-1.7.8.tar.gz
