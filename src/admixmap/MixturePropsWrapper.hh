#ifndef MIXTUREPROPSWRAPPER_H
#define MIXTUREPROPSWRAPPER_H
#include "utils/GSLErrorHandler.h"

/// wrapper class to provide an interface for the HMM to the mixture proportions
class MixturePropsWrapper{
public:
  MixturePropsWrapper(){
    theta = 0;
    stride = 0;
  }

  MixturePropsWrapper(const double* const p, unsigned K = 0){
    assign(p, K);
  }

  void assign(const double* const p, unsigned K = 0){
    theta = p;
    stride = K;
  }

  ~MixturePropsWrapper(){
    theta = 0;
    stride = 0;
  }

  double get(unsigned locus, unsigned state)const{
    if(!theta) 
      //      throw ("Error in MixturePropsWrapper: null pointer");
      GSLErrorHandler("Null pointer", __FILE__, __LINE__, -1);

    return theta[locus*stride + state];
  }

  double operator()(unsigned locus, unsigned state)const{
    return get(locus, state);
  }

  bool isNull()const{
    return (theta == NULL);
  }

  MixturePropsWrapper operator+(unsigned t)const{
    MixturePropsWrapper MPW(theta+t , stride);
    return MPW;
  }

private:
  const double* theta;
  unsigned stride;

};

#endif
