//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file Chromosome.cc
/// Implementation of the Chromosome class.
//=============================================================================

#include "Chromosome.h"
#include <iostream>
#include "bclib/misc.h"
//#include <string>    //for throwing exceptions in f
//#include <exception> // for catching exceptions in f

#define PR(x) cerr << #x << " = " << x << endl;

using namespace std;

/// @parm inNumHiddenStates Not really the number of hidden states.
Chromosome::Chromosome(int n, int size, int start, int inNumHiddenStates, bool isx = false):
  NumberOfCompositeLoci(size), isX(isx), Number(n), _startLocus(start), NumHiddenStates(inNumHiddenStates)
 {
   
   if ( size < 1 ){
     size = 1;
   }
   Distances = new double[ NumberOfCompositeLoci ];
   
   CodedStates = new int[size];
   f = new double[2*size];
   f[0] = f[1] = 0.0;
   
   HMM = new HiddenMarkovModel( size, NumHiddenStates, NumHiddenStates*NumHiddenStates, f);
 }

Chromosome::~Chromosome()
{
  delete[] CodedStates;
  delete[] f;
  delete[] Distances; 
  delete HMM;
}


// ******** Chromosome information **********************************
bool Chromosome::isXChromosome()const
    {
    return isX;
    }


/// Sets label for this chromosome, usually "1", "2", etc or "X"
void Chromosome::SetLabel( std::string label )
    {
    _Label = label;
    }


const string Chromosome::GetLabel( )const
    {
    return _Label;
    }


void Chromosome::SetDistance( unsigned int locus, double distance )
    {
    gp_assert_gt( locus, 0			    );
    gp_assert_lt( locus, GetNumberOfCompositeLoci() );

    Distances[ locus ] = distance;
    }


int Chromosome::GetLocus(int num)const{
  return _startLocus + num;
}



// ****************** Setting of locus correlation, f *************************

///returns locus correlation at a given locus, given a value of sum-intensities
//TODO: possibly make inline
double Chromosome::LocusCorrelation(unsigned locus, double drho){
  //  try{
    if(isX)//sumintensities on Xchrm is set to half autosomal value
      return bclib::eh_exp(-GetDistance( locus ) * drho*0.5);
    else
      return bclib::eh_exp(-GetDistance( locus ) * drho);
//  }
//   catch(std::exception e){
//     std::string err = "Error encountered while setting locus f in Chromosome: ";
//     err.append(e.what());
//     throw err;
//   }
}

/**
   Sets locus correlation, f.
   rho should be either
   (1) vector of length 1 containing a global (global rho model) or 
   individual-level (non globalrho and assortative mating) sumintensities
   (2) vector of length 2 containing gamete-specific sumintensities for a single individual
*/
void Chromosome::SetLocusCorrelation( const genepi::RhoType & vrho, bool RandomMating ){
  for( unsigned int j = 1; j < NumberOfCompositeLoci; j++ ){
    //first gamete
    f[2*j] = LocusCorrelation(j, vrho[0]);
    //second gamete
    if(RandomMating)//value for second gamete is different
      f[2*j + 1] = LocusCorrelation(j, vrho[1]);
    else//second gamete is the same
      //SetLocusCorrelelation(j, vrho[0]);
      f[2*j + 1] = f[2*j];
  }
}
void Chromosome::SetGlobalLocusCorrelation(double drho){
  for(unsigned j = 1; j < NumberOfCompositeLoci; j++ ){
    f[2*j    ] = LocusCorrelation(j, drho);
    f[2*j + 1] = f[2*j];
  }
}

// ********** Interface to HMM ****************************************
/**
   sets conditional probabilities of ancestry at locus j.
   One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
   i.e. (i,2) = p_{ii}

        (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}

	(i,0) = 1.0 - (i,1) - (i,2)
	where p's are probs in StateProbs from HMM
*/
void Chromosome::getHiddenStateCopyNumberProbs(std::vector<std::vector<double> >& CopyNumberProbs, const bool isDiploid, int j) {
  const bclib::pvector<double>& probs = HMM->GetHiddenStateProbs(isDiploid, j);

  if(isDiploid){
    for( int k1 = 0; k1 < NumHiddenStates; k1++ ){
      // prob of 2 copies = diagonal element
      CopyNumberProbs[2][k1] = probs[(NumHiddenStates + 1) * k1];
      CopyNumberProbs[1][k1] = 0.0;
      for( int k2 = 0 ; k2 < NumHiddenStates; k2++ )
        CopyNumberProbs[1][k1] += probs[k1*NumHiddenStates + k2]
                                + probs[k2*NumHiddenStates + k1];
      CopyNumberProbs[1][k1] -= 2.0*CopyNumberProbs[2][k1];
      CopyNumberProbs[0][k1] = 1.0 - CopyNumberProbs[1][k1] - CopyNumberProbs[2][k1];
    }
  }
  else{//haploid case
    for( int k = 0; k < NumHiddenStates; k++ ){
      CopyNumberProbs[2][k] = 0.0;            // haploid: cannot have two copies
      CopyNumberProbs[1][k] = probs[k];       // probability of one copy
      CopyNumberProbs[0][k] = 1.0 - probs[k]; // probability of no copies
    }
  }
}

///samples jump indicators xi for this chromosome and updates SumLocusAncestry
void Chromosome::SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
				      int *SumLocusAncestry, std::vector<unsigned> &SumN, bool SampleArrivals)const {
  HMM->SampleJumpIndicators(LocusAncestry, gametes, SumLocusAncestry, SumN, SampleArrivals, _startLocus);
}
