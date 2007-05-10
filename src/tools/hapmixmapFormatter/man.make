EXEC_NAME=FPHD
CC =g++
CXXFLAGS =-O3
BAYESLIB_PATH =../../bcppcl
INCLUDES =-I$(BAYESLIB_PATH)/utils

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl
#GSL_PATH
#LIBS =-L$(GSL_PATH)

SUFFIXES =.cc.cpp

objects =FormatPhasedData.o GenotypeEncoding.o FPHDOptions.o HapMapLegend.o 

all: $(objects)
	$(CXX) $(CXXFLAGS) -o$(EXEC_NAME) $(objects) $(LIBS) -lbcppcl -lgsl

%.o: %.cpp#rule for compilation
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc#rule for compilation
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o $(EXEC_NAME)