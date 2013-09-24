

LIB=-lIL


CFLAGS=-fopenmp
CXXFLAGS=${CFLAGS}
CXX=g++

all: z1

obj/z1.o: src/z1.cpp
	${CXX} -c ${CXXFLAGS} src/z1.cpp -o obj/z1.o

obj/camera.o: src/Camera.cpp src/Camera.hpp
	${CXX} -c ${CXXFLAGS} src/Camera.cpp -o obj/camera.o

z1: obj/z1.o obj/camera.o
	${CXX}  ${CXXFLAGS} obj/z1.o obj/camera.o -o bin/z1 ${LIB}
