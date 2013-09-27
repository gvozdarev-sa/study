

LIB=-lIL
#-ljasper -lmng -lpng -ljpeg -ltiff -Iinclude -Llib


CFLAGS=-O2
CXXFLAGS=${CFLAGS}
OPENMP=-fopenmp
CXX=g++

#########################
TASK=z2
#########################
all: Z1

Z1: clean run


obj/image.o: src/Image.cpp src/Image.hpp
	${CXX} -c ${CXXFLAGS} src/Image.cpp -o obj/image.o

### with -fopenmp
obj/z.o: src/${TASK}.cpp
	${CXX} -c ${CXXFLAGS} ${OPENMP} src/${TASK}.cpp -o obj/z1.o

bin/z: obj/z.o obj/image.o
	${CXX}  ${CXXFLAGS} ${OPENMP} obj/z1.o obj/image.o -o bin/z ${LIB}

z: bin/z

### w/0 -fopenmp
obj/zw.o: src/${TASK}.cpp
	${CXX} -c ${CXXFLAGS} src/${TASK}.cpp -o obj/zw.o

bin/zw: obj/zw.o obj/image.o
	${CXX}  ${CXXFLAGS}  obj/zw.o obj/image.o -o bin/zw ${LIB}

zw: bin/zw
##############


run: zw z
	cd bin && ./z 1.jpg
	cd bin && rm 2.bmp
	cd bin && ./zw 1.jpg

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
