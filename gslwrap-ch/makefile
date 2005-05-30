# makefile for creating vector and matrix library
CC 			= gcc
CXX 			= g++

LIBS			= -lm -static -lgsl -lgslcblas

CXXFLAGS		= 
CPPFLAGS		= -W -Wall -O3 -march=i686
#CPPFLAGS		= -W -Wall -g -march=i686

all:			libmat.a

OBJS			= matrix.o vector.o matrix_d.o vector_d.o matrix_i.o vector_i.o rand_serial.o

libmat.a:	$(OBJS)
		ar rv libmat.a $(OBJS)

matrix:		
	perl generate

clean:
	- rm *.o
	- rm *.a

