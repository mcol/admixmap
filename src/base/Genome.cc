//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
// Portions copyright (C) 2009  David D. Favro
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
/// \file Genome.cc
/// Implementation of the Genome class.
//=============================================================================

#include "Genome.h"
#include "Chromosome.h"
#include "SimpleLocusArray.h"

#include <cerrno>
#include <cstring>	// strerror()
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "bclib/DataMatrix.h"
#include "bclib/StringConvertor.h"
#include "bclib/LogWriter.h"
#include "bclib/estr.h"


using namespace std;
using namespace genepi;



/// This just exists so that we don't need to include <stdexcept> in the header
/// to throw exceptions in inline methods.
void Genome::throwErr( const std::string & msg ) const
    {
    throw runtime_error( msg );
    }



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  NumberOfChromosomes = 0;
  LocusArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  XChromosomeIndex = 0; 
  TotalLoci = 0;
  Distances = 0;
  SizesOfChromosomes = 0;
  C = 0;
}

Genome::~Genome()
{
  delete[] LocusArray;
  delete[] SizesOfChromosomes;
  delete[] Distances; 
  if(C){
      for(unsigned i = 0; i < NumberOfChromosomes; i++){
	  delete C[i];
      }
      delete[] C;
  }
}

///gets contents of locusfile and creates CompositeLocus array and Chromosome array
void Genome::Initialise(const InputData* const data_, int populations,
                        bool /*hapmixmodelindicator*/, bclib::LogWriter& Log) {

  const SimpleLocusArray & simpleLoci = data_->getSimpleLoci();
  
  //determine number of composite loci
  NumberOfCompositeLoci = data_->getNumberOfCompositeLoci();
  TotalLoci = data_->getNumberOfSimpleLoci();
  
  //create array of CompositeLocus objects
  LocusArray = new CompositeLocus[ NumberOfCompositeLoci ];
  Distances = new double[ NumberOfCompositeLoci ];
  
  // Set number of alleles at each locus
  unsigned row = 0; // counts lines in locusfile
  vector<unsigned int> cstart; ///< The starting composite-locus-index of each chromosome
  vector<size_t>       sstart; ///< The starting simple-locus-index of each chromosome
 
  int cnum = -1; //cnum = number of chromosomes -1
  int lnum = 0;
  LocusTable.resize(NumberOfCompositeLoci);
  X_data = false;

  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    LocusTable[i].resize(2);
    
    //retrieve first row of this comp locus from locusfile
    const SimpleLocus & sLoc = simpleLoci[row];
    if ( ! sLoc.startsNewChromosome() )
      Distances[ i ] = sLoc.getDistance().inMorgans();

    if ( sLoc.startsNewChromosome() ) {
      //new chromosome, triggered by missing value or distance of >= threshold
      cnum++;
      lnum = 0; 
      cstart.push_back(i  ); // composite-locus index of first locus on new chromosome
      sstart.push_back(row); // simple-locus index of first locus on new chromosome
    } else{
      lnum++;//one more locus on chromosome
    }
    LocusTable[i][0] = cnum;//chromosome on which locus i is located
    LocusTable[i][1] = lnum;//number on chromosome cnum of locus i

    // Set number of alleles of first locus in comp locus
    LocusArray[i].AddLocus( sLoc.getNumAlleles(), sLoc.getName() );

    // Loop through lines in locusfile for current complocus
    const SimpleLocus * nextLoc;
    while ( (++row != simpleLoci.size()) &&
            (nextLoc = &(simpleLoci[row]))->isCompositeWithPrevious() ) {
      // Adds locus with number of alleles given as argument:
      LocusArray[i].AddLocus( nextLoc->getNumAlleles(), nextLoc->getName() );
    }

    if(LocusArray[i].GetNumberOfLoci()>8) Log << "WARNING: Composite locus with >8 loci\n";
    //TotalLoci += LocusArray[i].GetNumberOfLoci();

  }//end comp loci loop

  NumberOfChromosomes = cnum +1;
  gp_assert_eq( NumberOfChromosomes, simpleLoci.getNChromosomes() );

  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];//array to store lengths of the chromosomes
  cstart.push_back(NumberOfCompositeLoci);//add extra element for next line to work
  for(unsigned c = 0; c < NumberOfChromosomes; ++c) SizesOfChromosomes[c] = cstart[c+1] - cstart[c];

  // create Chromosome objects
  InitialiseChromosomes(cstart, sstart, populations, simpleLoci);

  //print length of genome, num loci, num chromosomes
  PrintSizes(Log, data_->getUnitOfDistanceAsString());
}

///Creates an array of pointers to Chromosome objects and sets their labels.
///Also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
///LengthOfXChrm 
void Genome::InitialiseChromosomes(const std::vector<unsigned>& cstart,
                                   const std::vector<size_t>& sstart,
                                   int populations,
                                   const SimpleLocusArray& sLoci) {

  simpleLoci = &sLoci;
  C = new Chromosome*[NumberOfChromosomes]; 
  //C is an array of chromosome pointers

  for(unsigned i = 0; i < NumberOfChromosomes; i++){//loop over chromsomes

    const size_t cStLocIdx = cstart[i];
    const size_t size = cstart[i+1] - cStLocIdx; //number of loci on chromosome i

    const size_t        sStLocIdx  = sstart[i];
    const SimpleLocus & startLocus = sLoci[ sStLocIdx ];
    const bool          isX        = startLocus.isXChrom() == CHR_IS_X;
    const string &      label      = startLocus.getChromLabel();

    if(isX){
      X_data = true;
      XChromosomeIndex = cStLocIdx;//index of first locus on X chromosome
    }

    CreateChromosome(i, size, isX, cStLocIdx, populations, label);
  }
}

/// Create a chromosome.
/// \param i       index of this chromosome
/// \param size    number of composite loci on this chromosome
/// \param isX     indicates if it the X chromosome
/// \param cstart  index of the first locus
/// \param pops    number of populations
/// \param label   label for the chromosome (usually "1", "2", etc or "X")
void Genome::CreateChromosome(unsigned i, unsigned size, bool isX,
                              unsigned cstart, int pops, const string& label) {
  // C[i] is a pointer to Chromosome
  C[i] = new Chromosome(i, size, cstart, pops, isX);
  C[i]->SetLabel(label);

  // loop over loci on chromosome
  for (size_t j = 1; j < size; ++j) {
    double distance = GetDistance(cstart + j);
    C[i]->SetDistance(j, distance);
    if (!isX) // length of genome does not include X chromosome
      LengthOfGenome += distance;
    else
      LengthOfXchrm += distance;
  }
}

///accesses the entire chromosome array
const Chromosome* const* Genome::getChromosomes()const{
  return C;
}
///accesses a chromosome
Chromosome *Genome::getChromosome(unsigned j){
  return C[j];
}

///accesses a composite locus
CompositeLocus* Genome::operator() ( int ElementNumber ) const
{
  if ( ElementNumber >= (int)NumberOfCompositeLoci ){
    cout << "WARNING: Genome::operator() Element Number";
    cout << ElementNumber << " > NumberOfCompositeLoci: " << NumberOfCompositeLoci << " accessed." << endl;
  }
  if ( ElementNumber < 0){
    cout << "WARNING: Genome::operator() ElementNumber";
    cout << ElementNumber << " < 0 accessed." << endl;
  }

  return &(LocusArray[ElementNumber]);
}

const CompositeLocus* Genome::GetLocus(int ElementNumber)const{
  return this->operator()(ElementNumber);
}



//-----------------------------------------------------------------------------
// operator[] (Locus array access)
//-----------------------------------------------------------------------------

CompositeLocus & Genome::operator[]( unsigned int locusIdx )
    {
    if ( locusIdx >= NumberOfCompositeLoci )
	throw invalid_argument( genepi::estr("locus idx ") + locusIdx +
			" exceeds number of loci " + NumberOfCompositeLoci );

    return LocusArray[ locusIdx ];
    }


/// Writes numbers of loci and chromosomes and length of genome to Log and screen.
/// unit is the unit of measurement of the distances in the locusfile (Morgans/centiMorgans) 
void Genome::PrintSizes(bclib::LogWriter &Log, const string& distanceUnit)const{
  
  Log.setDisplayMode(bclib::Quiet);
  Log << "\n" << (int)TotalLoci << " simple loci\n"
      << (int)NumberOfCompositeLoci << " compound loci; "
      << (int)NumberOfChromosomes << " chromosome"; if(NumberOfChromosomes > 1) Log << "s";
  Log << "\n";

  Log << "Effective length of autosomes under study: ";
  if(distanceUnit == "cM")Log << LengthOfGenome*100.0 ;
  else Log << LengthOfGenome;
  Log << distanceUnit << ".\n";

  if( isX_data() ){
    Log << "Effective length of X chromosome under study: ";
    if(distanceUnit == "cM")Log << LengthOfXchrm*100.0;
    else Log << LengthOfXchrm;
    Log << distanceUnit << ".\n";


   }
  Log << "\n";
}

//Accessors

/// returns the number of composite loci
unsigned int Genome::GetNumberOfCompositeLoci()const
{
  return NumberOfCompositeLoci;
}
///returns the number of loci in a given composite locus
int Genome::getNumberOfLoci(int j)const{
  return LocusArray[j].GetNumberOfLoci();
}
///returns the number of chromosomes
unsigned int Genome::GetNumberOfChromosomes()const{
  return NumberOfChromosomes;
}

///returns total number of simple loci
unsigned int Genome::GetTotalNumberOfLoci()const{
  return TotalLoci;
}
///returns int array of chromosome sizes
const unsigned int *Genome::GetSizesOfChromosomes()const{
  return SizesOfChromosomes;
}

#if ALLOW_UNSAFE_DISTANCE_ACCESS
    /// returns the vector of distances between loci
    const double *Genome::GetDistances()const
    {
      return( Distances );
    }
#endif

///returns distance between a given locus and the previous one
double Genome::GetDistance( int locus ) const
{
  #if AGGRESSIVE_RANGE_CHECK

      if ( (locus < 0) || (unsigned(locus) > NumberOfCompositeLoci) )
	  throw runtime_error( estr("Access of out-of-bounds distance[") +
				locus + "] (max " + NumberOfCompositeLoci + ')' );

      if ( (*simpleLoci).size() == (*simpleLoci).getNComposite() )
	  gp_assert( ! (*simpleLoci)[locus].startsNewChromosome() );

  #endif

  return( Distances[ locus ] );
}

///returns total number of states accross all comp loci
int Genome::GetNumberOfStates()const{
  int ret = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += LocusArray[i].GetNumberOfStates();
  }
  return ret;
}

/// Get a chromosome number by absolute locus number
unsigned Genome::GetChrNumOfLocus(unsigned locus) const {
  return LocusTable[locus][0];
}

/// Get a chromosome-relative locus number by absolute locus number
int Genome::getRelativeLocusNumber(int j)const{
  return LocusTable[j][1];
}

/** Returns a vector of length 2
 * With chromosome number on 1st position and locus index on 2nd.
 */
const vector<int> Genome::GetChrmAndLocus(int j)const{
  return LocusTable[j];
}

const vector<vector<int > > Genome::GetChrmAndLocus()const{
  return LocusTable;
}
/// For a given (composite) locus, returns the number of the chromosome it is on and the position on that chromosome.
void Genome::GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l){
  *c = LocusTable[locus][0];
  *l = LocusTable[locus][1];
}
///indicates whether there is an X chromosome
bool Genome::isX_data()const
{
   return X_data;
}
///returns length of Genome in Morgans
double Genome::GetLengthOfGenome()const
{
   return LengthOfGenome;
}
//returns length of X chromosome in Morgans
double Genome::GetLengthOfXchrm()const
{
   return LengthOfXchrm;
}
///indicates if a chromosome is an X chromosome
unsigned Genome::isXChromosome(unsigned j)const{
  return C[j]->isXChromosome();
}
bool Genome::isXLocus(unsigned j)const{
  return C[LocusTable[j][0]]->isXChromosome();
}
/// returns index of X chromosome
unsigned Genome::getFirstXLocus()const{
  if(X_data)return XChromosomeIndex;
  else return NumberOfCompositeLoci;
}



//-----------------------------------------------------------------------------
/// Set global locus correlation across all chromosomes, case of vector-valued rho
/// In case of global rho model (rho has length 1), sets f globally across loci.
/// <B>NB:</B> <I>if @a rho has length other than 1, it is <B>ignored</B></I>
//-----------------------------------------------------------------------------

void Genome::SetLocusCorrelation( const genepi::cvector<double> & rho )
    {
    if ( rho.size() == 1 )
	for ( unsigned int j = 0; j < NumberOfChromosomes; j++ )
	    C[j]->SetLocusCorrelation(rho, false/*<-no random-mating*/);
    }



///set global locus correlation across all chromosomes, case of global rho
void Genome::SetLocusCorrelation(double rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->SetGlobalLocusCorrelation(rho);
  }
}



/// Prints table of composite loci for R script to read
void Genome::PrintLocusTable( const char * filename, const SimpleLocusArray & sLoci ) const
    {

    // Compilers other than gcc rarely support nested functions, yet we can
    // nest a class definition with an inline method here.
    class non_gcc {
	public:
	static inline void print( ostream & os,
				  const CompositeLocus & compLoc,
				  const string & cLabel, double mapPos )
	    {
	    os << compLoc.GetLabel(0) << '\t'
	       << compLoc.GetNumberOfStates() << '\t'
	       << mapPos << '\t' << cLabel << endl;
	    }
	};


    ofstream os( filename );
    if ( os.bad() || (! os.is_open()) )
	throw runtime_error( string( "failed to open file \"") +
			     filename + "\" because: " + strerror(errno) );

    os.exceptions( ios_base::badbit );

    os << "LocusName\tNumHaps\tMapPosition(" << gduAsString(sLoci.getGDU())
       << ")\tChromosome\n";

    unsigned int locus        = 0; // counter for composite locus
    unsigned int simple_locus = 0; // need to count simple loci to step through vector of distances

    for ( unsigned int c = 0 ; c < NumberOfChromosomes ; ++c )
	{
	const string & label = C[c]->GetLabel();
	double mapPosition   = 0.0;

	// Print the first locus on chromosome:
	non_gcc::print( os, LocusArray[locus], label, mapPosition );

	simple_locus += LocusArray[locus].GetNumberOfLoci();
	++locus;

	// Step through rest of loci on chromosome:
	for ( unsigned j = 1; j < SizesOfChromosomes[c]; ++j )
	    {
	    // Increment map position by distance of first locus in complocus:
	    //mapPosition += sLoci[simple_locus].getDistance().inCentimorgans();
	    mapPosition += sLoci[simple_locus].getDistance().inMorgans();
	    non_gcc::print( os, LocusArray[locus], label, mapPosition );
	    simple_locus += LocusArray[locus].GetNumberOfLoci();
	    ++locus;
	    }
	}
    }
