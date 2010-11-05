//=============================================================================
//
// Copyright (C) 2007  David O'Donnell
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file GenotypeIterator.h
/// Definition of the GenotypeIterator class.
//=============================================================================

#ifndef GENOTYPE_ITERATOR_H
#define GENOTYPE_ITERATOR_H


#include "config.h" // USE_GENOTYPE_PARSER
#if USE_GENOTYPE_PARSER
    // Included just for GType definition:
    #include "GFileLexer.h"
#else
    #include <vector>
#endif


/** \addtogroup base
 * @{ */


enum ploidy { haploid, diploid };


/// This is not an iterator class in the standard sense of the word, but rather
/// a sort of abstract "array adapter" into a genome (of varying type)
/// object that can be <I>used</I> to iterate via the "subscripts" to
/// operator(), j and g.
class GenotypeIterator {
 public:
  virtual ~GenotypeIterator() {}
  //virtual std::vector<unsigned short>::const_iterator operator[](unsigned j)const = 0;
  virtual unsigned short operator()(unsigned j, unsigned g) const = 0;
  unsigned short get(unsigned j, unsigned g) const { return operator()(j,g); }
  bool isDiploid() const { return isdiploid; }

  GenotypeIterator( bool _isdiploid ) : isdiploid( _isdiploid ) {}
  GenotypeIterator() {}

protected:
  bool isdiploid;
};


// For admixmap: one of two varieties, depending on whether we use old-style or
// new-style input-file parser:
#if USE_GENOTYPE_PARSER
    class AdmixGenotypeIterator : public GenotypeIterator {
     private:
      const genepi::GenotypeArray & v;
     public:
      AdmixGenotypeIterator(  const genepi::GenotypeArray & in, ploidy p) :
	GenotypeIterator( p == diploid ) ,
	v		( in	       ) {}
      unsigned short operator() (unsigned j, unsigned g) const;
      const genepi::Genotype & operator[] ( size_t j ) const { return v[ j ]; }
    };
#else
    class AdmixGenotypeIterator : public GenotypeIterator {
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
      //  std::vector<unsigned short>::const_iterator operator[](unsigned j)const{
      //  return (it+j)->begin();
      //}
     private:
      std::vector<std::vector<unsigned short> >::const_iterator it;
      AdmixGenotypeIterator();
    };
#endif


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
 private:
  std::vector<unsigned short>::const_iterator it;
  HapMixGenotypeIterator();
};


/** @} */


#endif
