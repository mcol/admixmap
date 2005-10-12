# makefile for creating vector and matrix library
CC 			= gcc
CXX 			= g++

LIBS			= -lm -static -lgsl -lgslcblas

CXXFLAGS		= 
CPPFLAGS		= -W -Wall -O3 -march=i686
#CPPFLAGS		= -W -Wall -g -march=i686

all:			libmat.a

OBJS			= rand_serial.o

libmat.a:	$(OBJS)
		ar rv libmat.a $(OBJS)

clean:
	- rm *.o
	- rm *.a

