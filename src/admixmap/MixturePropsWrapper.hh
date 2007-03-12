#ifndef MIXTUREPROPSWRAPPER_H
#define MIXTUREPROPSWRAPPER_H

/// wrapper class to provide an interface for the HMM to the mixture proportions
class MixturePropsWrapper{
public:
  MixturePropsWrapper(){
    theta = 0;
    K = 0;
  }

  MixturePropsWrapper(const double* const p, unsigned NumStates = 0){
    theta = p;
    K = NumStates;
  }

  ~MixturePropsWrapper(){
    theta = 0;
    K = 0;
  }

  double get(unsigned locus, unsigned state)const{
    return theta[locus*K + state];
  }

  double operator()(unsigned locus, unsigned state)const{
    return get(locus, state);
  }

  bool isNull()const{
    return (theta == NULL);
  }

  MixturePropsWrapper operator+(unsigned t)const{
    MixturePropsWrapper MPW(theta+t , K);
    return MPW;
  }

private:
  const double* theta;
  unsigned K;

};

#endif
