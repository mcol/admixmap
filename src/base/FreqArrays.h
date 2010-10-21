// *-*-C++-*-*
//=============================================================================
/// \file FreqArrays.h
/// Definition of the FreqArrays class.
//=============================================================================

#ifndef FREQ_ARRAYS_H
#define FREQ_ARRAYS_H


/** \addtogroup base
 * @{ */


/**
 * struct to hold allelecounts in either a 1d (where Number of alleles is fixed) or 2d array. 
 * usage, 1d array:
 * array_of_allelecounts AlleleCounts;
 * AlleleCounts.stride = K*2;
 * AlleleCounts.array = new int[ NumberOfLoci*K*2];
 *
 * usage, 2darray:
 * AlleleCounts.array = new int*[NumberOfLoci ];
 * for(unsigned i = 0; i < NumberOfLoci; ++i) AlleleCounts.array[i] = new int[K*NumberOfStates[i]];
 *
 * AlleleCounts[i]; //accesses counts for ith locus
 *
 * AlleleCounts[i][k*2 +a]; //accesses count of ath allele in kth pop at ith locus
 */

struct array_of_allelecounts {
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
};


/**
   struct to hold allelefreqs in either a 1d (where Number of alleles is fixed) or 2d array. 
   See array_of_allelecounts for details.
*/
class FreqArray
{

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



/** @} */


#endif
