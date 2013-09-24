

LIB=-lIL


CFLAGS=-fopenmp
CXXFLAGS=${CFLAGS}
CXX=g++

all: z1

obj/z1.o: src/z1.cpp
	${CXX} -c ${CXXFLAGS} src/z1.cpp -o obj/z1.o

obj/image.o: src/Image.cpp src/Image.hpp
	${CXX} -c ${CXXFLAGS} src/Image.cpp -o obj/image.o

z1: obj/z1.o obj/image.o
	${CXX}  ${CXXFLAGS} obj/z1.o obj/image.o -o bin/z1 ${LIB}
