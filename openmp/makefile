

LIB=-lIL
#-ljasper -lmng -lpng -ljpeg -ltiff -Iinclude -Llib


CFLAGS=-O2
CXXFLAGS=${CFLAGS}
OPENMP=-fopenmp
CXX=g++

all: Z1

Z1: clean run


obj/image.o: src/Image.cpp src/Image.hpp
	${CXX} -c ${CXXFLAGS} src/Image.cpp -o obj/image.o

### with -fopenmp
obj/z1.o: src/z1.cpp
	${CXX} -c ${CXXFLAGS} ${OPENMP} src/z1.cpp -o obj/z1.o

bin/z1: obj/z1.o obj/image.o
	${CXX}  ${CXXFLAGS} ${OPENMP} obj/z1.o obj/image.o -o bin/z1 ${LIB}

z1: bin/z1

### w/0 -fopenmp
obj/z1w.o: src/z1.cpp
	${CXX} -c ${CXXFLAGS} src/z1.cpp -o obj/z1w.o

bin/z1w: obj/z1w.o obj/image.o
	${CXX}  ${CXXFLAGS}  obj/z1w.o obj/image.o -o bin/z1w ${LIB}

z1w: bin/z1w
##############


run: z1w z1
	echo "Run with openmp\n============="
	cd bin && ./z1 1.jpg
	cd bin && rm 2.bmp
	echo "Run without openmp\n============"
	cd bin && ./z1w 1.jpg

clean: dirs
	rm -rf bin/*
	rm -rf obj/*
	cp src/1.jpg bin

dirs:
	mkdir -p bin
	mkdir -p obj
	mkdir -p lib
	mkdir -p include

devil: dirs
	wget http://downloads.sourceforge.net/openil/DevIL-1.7.8.tar.gz
	tar -xvf DevIL-1.7.8.tar.gz
	cd  devil-1.7.8 && ./configure
	cd  devil-1.7.8 && make
	cd  devil-1.7.8 && make install DESTDIR=`pwd`/..
	cp -rf usr/local/include/* include
	cp -rf usr/local/lib/* lib
	rm -rf devil-1.7.8
	rm -rf usr
	rm -f  lib/*.so*
	rm DevIL-1.7.8.tar.gz
