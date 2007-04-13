// *-*-C++-*-*
#ifndef HAPPAIR_H
#define HAPPAIR_H

#include <ostream>

///struct to hold a pair of haplotypes, coded as integers
class hapPair
 {
 public:
   int haps[2];

   hapPair(){
     haps[0] = haps[1] = -1;
   }
   hapPair(int a, int b){
     haps[0] = a;
     haps[1] = b;
   }

   friend 
   ///for printing a happair
   std::ostream& operator<<(std::ostream& os, const hapPair &h){
     os << h.haps[0] << " " << h.haps[1];
     return os;
   }

}; 


#endif
