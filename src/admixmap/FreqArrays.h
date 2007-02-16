// *-*-C++-*-*

#ifndef FREQ_ARRAYS_H
#define FREQ_ARRAYS_H

#include "config.h"
#include "utils/ColumnIter.h"

#ifndef PARALLEL
#define ARRAY2D
#endif

/**
   struct to hold allelecounts in either a 1d (where Number of alleles is fixed) or 2d array. 
   usage, 1d array:
   array_of_allelecounts AlleleCounts;
   AlleleCounts.stride = K*2;
   AlleleCounts.array = new int[ NumberOfLoci*K*2];
   
   usage, 2darray:
   AlleleCounts.array = new int*[NumberOfLoci ];
   for(unsigned i = 0; i < NumberOfLoci; ++i) AlleleCounts.array[i] = new int[K*NumberOfStates[i]];
   
   AlleleCounts[i]; //accesses counts for ith locus
   
   AlleleCounts[i][k*2 +a]; //accesses count of ath allele in kth pop at ith locus
   */
#ifdef ARRAY2D
typedef struct{

  int **array;

  int* operator[](unsigned i){//for reading or writing element
    return array[i];
  };
  const int* operator[](unsigned i)const{//for read-only
    return array[i];
  };
  void dealloc(int L){
    if(array){
      for(int i = L-1; i >=0 ; --i)
	if(array[i]){
	  delete[] array[i];
	  array[i] = 0;
	}
      delete[]array;
      array = 0;
    }
  };
}array_of_allelecounts;

/**
   struct to hold allelefreqs in either a 1d (where Number of alleles is fixed) or 2d array. 
   See array_of_allelecounts for details.
*/
class FreqArray {

public:
  double **array;

  FreqArray(){
   array = 0;
  };
  FreqArray(double** a){
   array = a;
  }

  FreqArray operator+(unsigned i){
   FreqArray A(array+i);
   return A;
  }
  double* operator[](unsigned i){//for reading or writing element
    return array[i];
  };
  const double* operator[](unsigned i)const{//for read-only
    return array[i];
  };
  void dealloc(int L){
    if(array){
      for(int i = L-1; i >=0 ; --i)
	if(array[i]){
	  delete[] array[i];
	  array[i] = 0;
	}
      delete[] array;
      array = 0;
    }
  };
};

#else
typedef struct{
  int* array;
  unsigned stride;

  int* operator[](unsigned i){
    return array + i*stride;
  };
  const int* operator[](unsigned i)const{
    return array + i*stride;
  };
  void dealloc(int ){
    if(array){
      delete[] array;
      array = 0;
    }
  };
}array_of_allelecounts;

class FreqArray {
public:
  double* array;
  unsigned stride;

  FreqArray(){
   array = 0;
   stride = 0;
  }

  FreqArray(double* a, unsigned s){
   array = a;
   stride = s;
  }

   ~FreqArray(){
   dealloc(0);
  }

  FreqArray operator+(unsigned i){
   FreqArray A(array+i, stride);
   return A;
  }
  double* operator[](unsigned i){
    return array + i*stride;
  };
  const double* operator[](unsigned i)const{
    return array + i*stride;
  };
  void dealloc(int ){
    if(array){
      delete[] array;
      array = 0;
    }
  };
};
#endif


#endif



