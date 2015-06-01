MATRICES_LSTDFLG = -lstdc++ -llapack -llapack_atlas -lcblas -lm -lgsl -lgslcblas
MATRICES_INCLUDE = -I/usr/include/ -I/usr/include/atlas/
MATRICES_LIB = -L/usr/lib/ -L/usr/lib/atlas-base/
MATRICES_OBJS = cov_test

all:	${MATRICES_OBJS}
	rm -f *.o

matrices.o: matrices.cpp matrices.h
	gcc -g -c matrices.cpp -o matrices.o ${MATRICES_INCLUDE}

hw7p1.o: hw7p1.cpp hw7p1.h
	gcc -g -c hw7p1.cpp -o hw7p1.o

main.o: main.cpp matrices.h hw7p1.cpp
	gcc -g -c main.cpp -o main.o ${MATRICES_INCLUDE}

cov_test: main.o matrices.o hw7p1.o
	gcc main.o matrices.o -o matrices ${MATRICES_LIB} ${MATRICES_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MATRICES_OBJS}
