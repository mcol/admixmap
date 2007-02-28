// *-*-C++-*-*
#ifndef GENOTYPE_ITERATOR_H
#define GENOTYPE_ITERATOR_H

enum ploidy {haploid, diploid};

class GenotypeIterator{
 public:
  virtual ~GenotypeIterator(){};
  //virtual std::vector<unsigned short>::const_iterator operator[](unsigned j)const = 0;
  virtual unsigned short operator()(unsigned j, unsigned g)const = 0;
  virtual unsigned short get(unsigned j, unsigned g)const = 0;
  bool isDiploid()const{return isdiploid;};

protected:
  bool isdiploid;

};

class AdmixGenotypeIterator : public GenotypeIterator{
 public:
  AdmixGenotypeIterator(  const std::vector<std::vector<unsigned short> >::const_iterator in, ploidy p){
    isdiploid = (p == diploid);
    it = in;
  }
  AdmixGenotypeIterator(  const std::vector<std::vector<unsigned short> >& in, ploidy p){
    isdiploid = (p == diploid);
    it = in.begin();
  }
  unsigned short operator()(unsigned j, unsigned g)const{
    if(g !=0 && g!= 1)throw ("Bad call to GenotypeOperator::operator()");
    return  ( *(it+j) )[g];
  }
  unsigned short get(unsigned j, unsigned g)const{
    return this->operator()(j,g);
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
    unsigned short a = 0;
    if(isdiploid){
      switch(*(it+j)){
      case 0: {
        a = 0;
        break;
      }
      case 1:{
        a = 1;
        break;
      }
      case 2:{
        a = 2;
        break;
      }
      case 3: {
        if(g==0) a = 1;
        else a = 2;
        break;
      }
      }
    }
    else{//haploid
      a = *(it+j);
    }
    return a;
  }
  unsigned short get(unsigned j, unsigned g)const{
    return this->operator()(j,g);
  }
 private:
  std::vector<unsigned short>::const_iterator it;
  HapMixGenotypeIterator();
};

#endif
