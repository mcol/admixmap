EXEC_NAME=FPHD
CC =g++
CXXFLAGS =-O3
BAYESLIB_PATH =../../bcppcl
INCLUDES =-I$(BAYESLIB_PATH)/utils#for StringSplitter

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to lik against entire gsl and cblas
#GSL_PATH
#LIBS =-L$(GSL_PATH)

SUFFIXES =.cc.cpp

objects =FormatPhasedData.o GenotypeEncoding.o FPHDOptions.o HapMapLegend.o $(BAYESLIB_PATH)/utils/StringSplitter.o

all: $(objects)
	$(CXX) $(CXXFLAGS) -o$(EXEC_NAME) $(objects) $(LIBS) -lgsl -lgslcblas

%.o: %.cpp#rule for compilation
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc#rule for compilation
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o $(EXEC_NAME)