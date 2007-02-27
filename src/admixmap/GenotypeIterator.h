// *-*-C++-*-*
#ifndef GENOTYPE_ITERATOR_H
#define GENOTYPE_ITERATOR_H

enum ploidy {haploid, diploid};

class GenotypeIterator{
 public:
  //virtual std::vector<unsigned short>::const_iterator operator[](unsigned j)const = 0;
  virtual unsigned short operator()(unsigned j, unsigned g)const = 0;
};

class AdmixGenotypeIterator : public GenotypeIterator{
 public:
  AdmixGenotypeIterator(  const std::vector<std::vector<unsigned short> >::const_iterator in){
    it = in;
  }
  AdmixGenotypeIterator(  const std::vector<std::vector<unsigned short> >& in){
    it = in.begin();
  }
  unsigned short operator()(unsigned j, unsigned g)const{
    if(g !=0 && g!= 1)throw ("Bad call to GenotypeOperator::operator()");
    return  ( *(it+j) )[g];
  }

  //  std::vector<unsigned short>::const_iterator operator[](unsigned j)const{
  //  return (it+j)->begin();
  //}
 private:
  std::vector<std::vector<unsigned short> >::const_iterator it;
  AdmixGenotypeIterator();

};

class HapMixGenotypeIterator : public GenotypeIterator{

 public:
  HapMixGenotypeIterator(const std::vector<unsigned short>::const_iterator in, ploidy p){
    it = in;
    isdiploid = (p == diploid);
  }
  unsigned short operator()(unsigned j, unsigned g = 0)const{
    if( (!isdiploid && g!=0) || (isdiploid && g!=1 && g!= 0) )throw ("Bad call to GenotypeOperator::operator()");
    if(isdiploid){
      switch(*(it+j)){
        case 0: return 0;
        case 1: return 1;
        case 2: return 2;
        case 3: {
	  if(g==0) return 1;
	  else return 2;
        }
      }
    }
    else{//haploid
      return *(it+j);
    }

  }

 private:
  std::vector<unsigned short>::const_iterator it;
  bool isdiploid;
  HapMixGenotypeIterator();
};

#endif
