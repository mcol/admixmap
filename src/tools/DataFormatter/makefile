EXEC_NAME=DataFormatter
CC =g++
CXXFLAGS =-O3
BAYESLIB_PATH =../../bayeslib
INCLUDES =-I$(BAYESLIB_PATH)/utils#for StringSplitter
#GSL_PATH
#LIBS =-L$(GSL_PATH)

SUFFIXES =.cc.cpp

objects =FormatPhasedData.o GenotypeEncoder.o Legend.o $(BAYESLIB_PATH)/utils/StringSplitter.o

all: $(objects)
	$(CXX) $(CXXFLAGS) -o$(EXEC_NAME) $(objects) $(LIBS) -lgsl -lgslcblas

%.o: %.cpp#rule for compilation
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cc#rule for compilation
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o $(EXEC_NAME)