//=============================================================================
//
// Copyright (C) 2007  David O'Donnell and Paul McKeigue
// Portions Copyright (C) 2009  David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 3 as published by the Free
// Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.	 If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file AdmixmapGenotypeConverter.cc
/// Implementation of the
/// genepi::convert(const Organism&org,const Genome&loci,std::vector<genotype> &genotypes,bool**missing)
/// function.
//=============================================================================

#include "AdmixmapGenotypeConverter.h"


#include <string>
#include <stdexcept>
#include <iostream>	// cerr

#include "bclib/estr.h"
#include "Genotype.h"


using namespace std;


/// Should we consider male with diploid genotype-data on X chromosome fatal error?
#define MALE_DIPLOID_X_FATAL	 0

/// Should male with diploid genotype-data on X chromosome issue warning if the
/// two gametes' alleles are the same?
#define MALE_DIPLOID_X_WARN_HOMO 0

/// Male diploid genotypes on X chromosome should be regarded as missing or haploid?
#define MALE_DIPLOID_X_HOMO_MISSING	0
#define MALE_DIPLOID_X_HETERO_MISSING	0

/// Should an individual with no genotypes be a fatal error?
#define NO_GENOTYPES_FATAL	 0


/// Not sure why we would, but we could count the number of loci on the X chromosome
#define COUNT_X_LOCI		 0


// If WAIT_DONT_THROW is turned on, when we encounter an error, just print a
// diagnostic to cerr, then throw one error after all processing is finished.
// If not, throw immediately.
#define WAIT_DONT_THROW		1



namespace genepi { // ----



typedef unsigned long Count;



static void thrGenErr( const Organism & org, const string & details, bool & haveErr )
    {
    const estr & msg =
	org.inLineDesc() + ": genotype error in organism " + org.idDesc() + ": " + details + '.';

    #if WAIT_DONT_THROW
	cerr << msg << '\n';
	haveErr = true;
    #else
	throw runtime_error( msg );
    #endif
    }



//-----------------------------------------------------------------------------
/// Validate the genotype input file for some obvious errors
//-----------------------------------------------------------------------------

void checkGenotypes( const Organism & org,
			bool haveObserved, Count numhaploid, Count numdiploid,
			Count numhaploidX, Count numdiploidX, bool & haveErr )
    {

    // Check for no observed genotypes (warning only)
    if ( ! haveObserved )
	{
	cerr << org.inLineDesc() << ": WARNING: individual " << org.idDesc()
		<< " has no observed genotypes." << endl;
	#if NO_GENOTYPES_FATAL
	    exit(1);
	#endif
	}

    // Check for male with diploid X data
    if ( (numhaploidX != 0) || (numdiploidX != 0) ) // Some X genotypes present
	{
	#if 0 // Apparently not checking for this
	    if ( (numdiploidX != 0) && !org.isFemale() )
		thrGenErr( ind, "males cannot have diploid X-chromosome genotypes" );
	#endif

	if ( (numhaploidX != 0) && org.isFemale() )
	    {

	    // Check for female with haploid and diploid X data
	    if ( numdiploidX != 0 )
		thrGenErr( org, "females should have diploid X-chromosome genotypes", haveErr );

	    // Check for phased X data but unphased autosomal genotypes
	    if ( numdiploid != 0 )
		thrGenErr( org, "female with diploid autosomes and haploid X-Chromosome.", haveErr );

	    }
	}

    else // only autosomes
	{
	// Check for mixed haploid/diploid data:
	if ( (numhaploid != 0) && (numdiploid != 0) )
	    thrGenErr( org, "both haploid and diploid genotypes and no X chromosome", haveErr );
	}

    }



//-----------------------------------------------------------------------------
/// Gets genotypes in admixmap model (hapmix genotypes are coded differently)
//-----------------------------------------------------------------------------

void convert( const Organism &		org	  ,
	      const Genome &		loci	  ,
	      vector<GenotypeArray> &	genotypes ,
	      bool**			missing	  )
    {

    SLocIdxType	 sLocIdx      = 0; // Simple locus counter
    unsigned int compLocusIdx = 0;
    Count numHaploid	= 0;
    Count numDiploid	= 0;
    Count numHaploidX	= 0;
    Count numDiploidX	= 0;
    #if COUNT_X_LOCI
	Count numXloci	= 0;
    #endif
    bool haveObserved = false;


    bool haveErr = false;

    unsigned int numberOfChromosomes = loci.GetNumberOfChromosomes();

    // Loop over each chromosome (c counter):
    for ( unsigned c = 0; c < numberOfChromosomes; ++c ) {

      const Chromosome & chromosome = loci.getChromosomeRef( c );
      const bool	 isXchrm    = chromosome.isXChromosome();

      // Within each chromosome, loop over the <<something>>
      const unsigned int cSize = loci.GetSizeOfChromosome( c );
      for ( unsigned int j = 0; j < cSize; ++j ) {

	GenotypeArray G;

	// Loop over composite loci to store genotypes in <<something>>

	const CompositeLocus &	compLocus = loci[ compLocusIdx ];
	const int		numLoci	  = compLocus.GetNumberOfLoci();
	bool			isMissing = true;

	#if COUNT_X_LOCI
	    if ( isXchrm )
		numXloci += numloci;
	#endif

	for ( int locus = 0; locus < numLoci; locus++ ) {

	  const GenotypeParser::GType & g = org.getGType( sLocIdx );

	  if ( isXchrm ) {
	    if ( g.isHaploid() ) {
	      ++numHaploidX;
	    }
	    else { // diploid X genotype
	      if ( ! org.isFemale() ) { // males cannot have diploid X genotypes
		// NOTE: allowing this for backward compatibility, for now
		// instead remove second element
		const bool heterozygous = g.hasTwoVals() && (g.getVal1() != g.getVal2());
		estr msg("at locus #");
		msg << sLocIdx << " ("
		    << org.getInFile().getSLoci()[sLocIdx].getName()
		    << "): male has diploid " << (heterozygous ? "heterozygous " : "")
		    << "X-chromosome genotype " << g.desc();
		#if MALE_DIPLOID_X_FATAL
		    thrGenErr( ind, msg, haveErr );
		#else

		    if ( heterozygous )
			#if MALE_DIPLOID_X_HETERO_MISSING
			    g.forceMissing();
			#else
			    g.forceHaploid();
			#endif
		    else
			#if MALE_DIPLOID_X_HOMO_MISSING
			    g.forceMissing();
			#else
			    g.forceHaploid();
			#endif

		    #if ! MALE_DIPLOID_X_WARN_HOMO
		      if ( heterozygous )
		    #endif
			{
			cerr << org.inLineDesc() << ": WARNING: individual " << org.idDesc() <<
			    ": " << msg
			#if MALE_DIPLOID_X_MISSING
				<< ": forcing to haploid\n";
			#else
				<< ": considering to be missing\n";
			#endif
			}

		#endif
	      }
	      ++numDiploidX;
	    }
	  }
	  else{//autosome
	    if ( g.isHaploid() )
		++numHaploid;
	    else ++numDiploid;
	  }
	  sLocIdx++;
	  G.push_back(g);

	  // If even a single simple-locus is not missing, the whole
	  // compound-locus is not missing; i.e. the compound locus is missing
	  // if all of its simple-loci are missing:
	  if ( ! g.isMissing2() )
	    isMissing = false;
	  haveObserved |= (! g.isMissing2());

	} // End locus loop

	missing[c][j] = isMissing;

	genotypes.push_back( G );
	++compLocusIdx;
      }

    }

    checkGenotypes( org, haveObserved, numHaploid, numDiploid, numHaploidX, numDiploidX, haveErr );

    if ( haveErr )
	throw runtime_error( "Fatal errors in genotype file" );
  }



} // ---- end namespace genepi
