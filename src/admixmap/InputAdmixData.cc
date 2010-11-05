//=============================================================================
//
// Copyright (C) 2007  David O'Donnell and Paul McKeigue
// Portions Copyright (C) 2009  David D. Favro
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
/// \file InputAdmixData.cc
/// Implementation of the InputAdmixData class.
//=============================================================================

#include "InputAdmixData.h"

#include "bclib/DataReader.h"
#include "config.h" // USE_GENOTYPE_PARSER
#include "bclib/estr.h"
#include "bclib/LogWriter.h"

#include <cstdlib>
#include <cstring>  // strlen()
#include <iostream>

#include "AlleleFreqParser.h"
#include "AlleleArray.h"

#include "AdmixOptions.h"
#include "CodeTimer.h"
#include "PedigreeAddl.h"
#include "AncestryVector.h" // For set_parms() [initialization]

#if USE_GENOTYPE_PARSER
    #include "AdmixmapGenotypeConverter.h"
#endif

#if defined(_OPENMP)
    #include <omp.h>
    #include "AdmixIndividualCollection.h"	// For PARALLELIZE_EPROB_COMPS
#endif

using bclib::LogWriter;
using namespace genepi;


#define DEBUG_PED_CONSTRUCTION	0



InputAdmixData::InputAdmixData( AdmixOptions & options, LogWriter & log )
    {
    using bclib::DataReader;

    #if ! USE_GENOTYPE_PARSER
      genotypeLoader = new GenotypeLoader;
    #endif

    log.setDisplayMode( bclib::Quiet );

    // Read all input files.
    Pedigree::setOptions( options );

    // Read base data files
    ReadData( &options, log );

    DataReader::ReadData( options.getAlleleFreqFilename(), alleleFreqData_, log );
    DataReader::ReadData( options.getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_, log );
    DataReader::ReadData( options.getEtaPriorFilename(), etaPriorData_,etaPriorMatrix_, log, false ); // no header
    DataReader::ReadData( options.getReportedAncestryFilename(), reportedAncestryData_, reportedAncestryMatrix_, log );

    log << "\n";

    CheckData( &options, log );
    }



//-------------------------------------------------------------------------
// This exists, due to the circular dependency between AdmixOptions and
// InputAdmixData.
//-------------------------------------------------------------------------

void InputAdmixData::finishConstructing( const AdmixOptions & options )
    {

    if ( isPedFile() || options.getUsePedForInd() )
	{


	int n_peds_excl = 0;


	#if ! (defined(NEEDED_ONLY_FOR_CONJUGATE_UPDATE) && NEEDED_ONLY_FOR_CONJUGATE_UPDATE)
	  if ( ! options.getNoConjugateUpdate() )
	    throw std::runtime_error( "Pedigrees are not compatible with conjugate-update."
					"  To correct, enable the no-conjugate-update option." );
	#endif


	// Generate pedigrees: must know how many populations, the _real_ value of
	// which isn't known until after CheckData(), so we delay until here
	// (i.e. otherwise we could and would call it from the constructor).
	generatePedigrees( options );


	// If the option is turned on, exclude pedigrees with more than the
	// specified number of members.
	if ( options.getExcludePedsOver() != 0 )
	    for ( size_t pIdx = getPeds().size() ; pIdx-- != 0 ; )
		if ( getPed(pIdx).getNMembers() > size_t(options.getExcludePedsOver()) )
		    {
		    cerr << "Excluding pedigree " << getPed(pIdx).getId()
			<< " because size=" << getPed(pIdx).getNMembers()
			<< " exceeds limit=" << options.getExcludePedsOver() << '\n';
		    getPeds().erase( pIdx );
		    }


	// AdmixedPedigree::InitialiseAdmixedStuff() needs to know the value of
	// options._admixed, which isn't known until after
	// InputAdmixOptions::checkOptions(), which isn't called [in admixmap.cc's
	// main()] until after InputAdmixData is finished constructing because it
	// needs the number of individuals passed into it.  This causes a nasty
	// circular dependency, and means that this must be delayed until after
	// InputAdmixOptions::checkOptions(), and can't be called from the
	// constructor of either this or the parent class.

	// Construct the new-format allele-frequencies matrix: for the moment, we'll
	// re-parse the input file (which was already read into one of
	// AdmixInputData's various matrices) rather than construct it from the
	// already-read data.

	AlleleProbVect		   alProbVect;
	AlleleFreqParser::PopArray populations;

	const PopIdx K = options.getPopulations();

	const char * const aFileName = options.getPriorAlleleFreqFilename();
	if ( (aFileName == 0) || (aFileName[0] == '\0') )
	    throw std::runtime_error( "at this time, pedigrees are only supported"
					" with a prior-allele-frequency file" );
	AlleleFreqParser::parse( aFileName, getSimpleLoci(), populations, alProbVect );
	gp_assert_eq( K, populations.size() );

	#if DEBUG_PED_CONSTRUCTION

	    fprintf( stderr, "\n\"Finishing\" intialize-peds, %zu pedigrees, %zu loci\n",
		getPeds().size(), getSimpleLoci().size() );

	    #if 1 // Turn on to print allele-frequency table
		cerr << "\n\n==== Allele-Probabilities ====\n\n";
		for ( SLocIdxType sLocIdx = 0 ; sLocIdx < getSimpleLoci().size() ; ++sLocIdx )
		    alProbVect[sLocIdx].print( cerr << '\n', populations );
		cerr << "\n\n";
	    #endif

	#endif


	// Big hack, we need the concept of 'context'. AncestryVector needs to
	// know some limits that are specific to the dataset.  This should
	// probably be replaced by something like an AncestryVectorDomain.  This
	// needs to be initialised prior to generating the hidden-state-space,
	// but after the pedigrees are constructed.

	Pedigree::FGameteIdx maxFounderGametes = 0;
	for ( int idx = 0 ; idx < int(getPeds().size()) ; ++idx )
	    {
	    // We assume that the number of founder-gametes on a non-X
	    // chromosome must be >= to that on an X chromosome; since we are
	    // only looking for the maxmimum, we only check that value.
	    const Pedigree::FGameteIdx nfg = getPed(idx).getNFounderGametes(CHR_IS_NOT_X);
	    if ( nfg > maxFounderGametes )
		maxFounderGametes = nfg;
	    }

	gp_assert( (getPeds().size() == 0) || (maxFounderGametes != 0) );
	gp_assert_le( maxFounderGametes, AV_MAX_FOUNDER_GAMETES );

	if ( options.getDisplayLevel() >= 3 )
	    cout << "Maximum number of founder-gametes: " << maxFounderGametes << '\n';

	AncestryVector::set_parms( options.getPopulations(), maxFounderGametes );


	genepi::CodeTimer ct;
	cerr << ct.local_started() << " computing hidden-state spaces and emission probabilities...\n";

	#if defined(_OPENMP) && PARALLELIZE_EPROB_COMPS
	  #pragma omp parallel for default(shared) EPROB_LOOP_OMP_SCHED if(options.getUsePedForInd())
	#endif
	  for ( int idx = 0 ; idx < int(getPeds().size()) ; ++idx )
	    {

	    Pedigree & ped = getPed( idx );

	    #if DEBUG_PED_CONSTRUCTION
		fprintf( stderr, "	ped %zu (%s) at %p: %zu members, %zu founders, %zu non-founders\n", idx,
		    ped.getId().c_str(), &ped, ped.getNMembers(), ped.getNFounders(), ped.getNNonFndrs() );
		fflush( stderr );
	    #endif

	    ped.genPossibleStatesInternal( K, alProbVect );

	    }


	// This can't be done in the parallelized loop above because
	// Mendelian-error-exclusions we mess with the iterator here.
	for ( int idx = 0 ; idx < int(getPeds().size()) ; ++idx )
	    {

	    Pedigree & ped = getPed( idx );

	    if ( options.getPrintPedSummary() )
		ped_sum_hss(cerr,ped) << '\n';

	    if ( options.getExcludeMendelError() )
		{

		// DDF: At this time, I am disqualifying all pedigrees which
		// contain mendelian errors/inconsistencies.  This is because
		// they will have emission probability of 0 for all hidden
		// states (i.e. hidden state space of size 0), because no
		// inheritance vector can ever be consistent with the observed
		// genotyped data.  If no correction is made, running the HMM
		// will generate errors because it will try to normalize
		// likelihoods that sum to 0.  We can add an exception to the
		// code that will deal with this, but haven't yet.

		if ( ped.getNMendelErrs() != 0 )
		    {
		    cerr << "Pedigree " << ped.getId() << " excluded due to "
			 << ped.getNMendelErrs() << " Mendelian errors!\n";
		    getPeds().erase( idx-- );
		    ++n_peds_excl;
		    continue;
		    }

		}


	    ped.InitialiseAdmixedStuff( options );

	    }


	if ( n_peds_excl != 0 )
	    cerr << "\nAfter excluding " << n_peds_excl
		 << " pedigrees due to Mendelian errors, "
		 << getPeds().size() << " pedigrees remain.\n";

	if ( getPeds().size() == 0 )
	    throw std::runtime_error( "No pedigree can be processed.\n" );

	cerr << ct.local_now() << " finished HSS/EP generation: " << ct.local_elapsed() << '\n';


	#if DEBUG_PED_CONSTRUCTION
	    fprintf( stderr, "---> Finished constructing InputAdmixData.\n" );
	#endif

	}

    }



//-----------------------------------------------------------------------------
// CheckData()
//-----------------------------------------------------------------------------

void InputAdmixData::CheckData(AdmixOptions *options, LogWriter &Log){

  #if ! USE_GENOTYPE_PARSER
    NumSimpleLoci = getNumberOfSimpleLoci();
    distanceUnit = DetermineUnitOfDistance();
    NumCompositeLoci = determineNumberOfCompositeLoci();
  #endif

  Log.setDisplayMode(bclib::Quiet);

  #if ! USE_GENOTYPE_PARSER
      bool badData = false;
      if(options->CheckData())
	badData = !checkLocusFile(options, Log);

      if(badData)
	exit(1);
  #endif

  SetLocusLabels();

  CheckAlleleFreqs(options, Log);
  ReadPopulationLabels(options);

  //detects regression model
  if(strlen( options->getOutcomeVarFilename() ) || strlen( options->getCoxOutcomeVarFilename() )){//if outcome specified
    const unsigned N = (genotypeLoader->getNumberOfIndividuals() - options->getTestOneIndivIndicator());
    if ( strlen( options->getOutcomeVarFilename() ) != 0 )
      CheckOutcomeVarFile( N, options, Log);
    if ( strlen( options->getCoxOutcomeVarFilename() ) != 0 ){
      OutcomeType.push_back( CoxData );
      if(options->CheckData())
	CheckCoxOutcomeVarFile( Log);
    }
    if ( strlen( options->getCovariatesFilename() ) != 0 )
      CheckCovariatesFile(N, options, Log);

    // append population labels to covariate labels
    if (!options->getTestForAdmixtureAssociation()) {
      for (vector<string>::const_iterator i = HiddenStateLabels.begin() + 1;
           i != HiddenStateLabels.end(); ++i) {
	// cout << "HiddenStateLabels " << *i << endl;
	CovariateLabels.push_back("slope." + *i); 
      }
    }
  }

  if ( strlen( options->getReportedAncestryFilename() ) != 0 )
    CheckRepAncestryFile(options->getPopulations(), Log);
}



//-----------------------------------------------------------------------------
// GetGenotype()
//-----------------------------------------------------------------------------

void InputAdmixData::GetGenotype(int i, const Genome &Loci, 
			   std::vector<genotype>* genotypes, bool **Missing) const {
  #if USE_GENOTYPE_PARSER
    convert( (*genotypeLoader)[i-1], Loci, *genotypes, Missing );
  #else
    genotypeLoader->GetGenotype(i, Loci, genotypes, Missing);
  #endif
}



//-----------------------------------------------------------------------------
// ReadPopulationLabels()
//-----------------------------------------------------------------------------

void InputAdmixData::ReadPopulationLabels(AdmixOptions *options){
  using bclib::DataReader;
  //  if(strlen(options->getAlleleFreqFilename()) || strlen(options->getPriorAlleleFreqFilename()) || strlen(options->getHistoricalAlleleFreqFilename())){
  if(strlen(options->getAlleleFreqFilename()))
    DataReader::ReadHeader(options->getAlleleFreqFilename(), HiddenStateLabels);
  else if ( strlen(options->getPriorAlleleFreqFilename()) != 0 )
    DataReader::ReadHeader(options->getPriorAlleleFreqFilename(), HiddenStateLabels);
  else if(strlen(options->getHistoricalAlleleFreqFilename()))
    DataReader::ReadHeader(options->getHistoricalAlleleFreqFilename(), HiddenStateLabels);
  else{
    //read labels from 'poplabels' option
    bclib::StringSplitter::Tokenize(options->getPopLabelString(), HiddenStateLabels, " ,");
    if(HiddenStateLabels.size() != (unsigned)options->getPopulations()){
      HiddenStateLabels.clear();

      //set default pop labels
      for( int j = 0; j < options->getPopulations(); j++ )
	HiddenStateLabels.push_back( genepi::estr("Pop") + (j+1) );
    }
  }
}
 

void InputAdmixData::CheckAlleleFreqs(AdmixOptions *options, LogWriter &Log){
  string freqtype = "";
  bool infile = false;//indicates whether either of the three allelefreq files are specified
  int nrows=0, expectednrows=0;
  int Populations = options->getPopulations();
  int NumberOfStates = 0;

  #if USE_GENOTYPE_PARSER

    int nStatesThisCompLoc = 0;
    for ( SimpleLocusArray::ConstIter slIter = getSimpleLoci().begin() ;
			slIter != getSimpleLoci().end() ; ++slIter )
	if ( slIter->isCompositeWithPrevious() )
	    nStatesThisCompLoc *= slIter->getNumAlleles();
	else
	    {
	    NumberOfStates += nStatesThisCompLoc;
	    nStatesThisCompLoc = slIter->getNumAlleles();
	    }
    NumberOfStates += nStatesThisCompLoc;

  #else // ! USE_GENOTYPE_PARSER:

    unsigned index = 0;
    for ( unsigned i = 0; i < NumCompositeLoci; ++i ) {
      int states = 1;
      do{
	states *= (int)locusMatrix_.get( index, 0 );
	index++;
	}
	while( index < locusMatrix_.nRows() && !locusMatrix_.isMissing(index, 1) && locusMatrix_.get( index, 1 ) == 0 );
      NumberOfStates += states;
    }

  #endif // ! USE_GENOTYPE_PARSER


  //fixed allele freqs
  if( strlen( options->getAlleleFreqFilename() )){
    freqtype = "";
    infile = true;
    nrows = alleleFreqData_.size()-1;
    expectednrows = NumberOfStates-getNumberOfCompositeLoci();
    Populations = alleleFreqData_[0].size() - 1;// -1 for ids in first col
    //getPopLabels(alleleFreqData_[0], Populations, PopulationLabels);
  }
  
  //Historic allelefreqs
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    freqtype = "historic";
    infile = true;
    nrows = historicalAlleleFreqData_.size();
    expectednrows = NumberOfStates+1;
    Populations = historicalAlleleFreqData_[0].size() - 1;
    //getPopLabels(historicalAlleleFreqData_[0], Populations, PopulationLabels);

  }
  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
    freqtype = "prior";
    infile = true;
    nrows = getPriorAlleleFreqData().size();
    expectednrows = NumberOfStates+1;
    Populations = getPriorAlleleFreqData()[0].size() - 1;
    //getPopLabels(getPriorAlleleFreqData()[0], Populations, PopulationLabels);
  }
  if(infile){
    if(nrows != expectednrows){
      Log << "Incorrect number of rows in " << freqtype << "allelefreqfile.\n" 
	  << "Expecting " << expectednrows << " rows, but there are " << nrows << " rows.\n";
      exit(1);
    }
    options->setPopulations(Populations);
  }
  else{//'populations' option
    if(Populations < 1){
      Log << "ERROR: populations = " << options->getPopulations() << "\n";
      exit(1);
    }

    //     for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    //       if(Loci->GetNumberOfStates(i) < 2){
    // 	Log << "ERROR: The number of alleles at a locus is " << Loci->GetNumberOfStates(i) << "\n";
    // 	exit(1);
    //       }
    //     }
  }

  genepi::Pedigree::setK( Populations ); // We need the concept of context
}

void InputAdmixData::CheckRepAncestryFile(int populations, LogWriter &Log)const{
  if( reportedAncestryMatrix_.nRows() != 2 * genotypeLoader->getNumberOfIndividuals() ){
    Log << "ERROR: " << "ReportedAncestry file has " << reportedAncestryMatrix_.nRows() << " rows\n"
	    "Genotypesfile has " << genotypeLoader->getNumberOfIndividuals() << " rows\n";
    exit(1);}
  if( (int)reportedAncestryMatrix_.nCols() != populations ){
    Log << "ERROR: " << "ReportedAncestry file has " << reportedAncestryMatrix_.nCols() << " cols\n"
	<< "AlleleFreq file has "<< populations << " cols\n";
    exit(1);
  }
}

const Matrix_s& InputAdmixData::getAlleleFreqData() const{
  return alleleFreqData_;
}

const Matrix_s& InputAdmixData::getHistoricalAlleleFreqData() const{
  return historicalAlleleFreqData_;
}

const Matrix_s& InputAdmixData::getEtaPriorData() const{
  return etaPriorData_;
}

const Matrix_s& InputAdmixData::getReportedAncestryData() const{
  return reportedAncestryData_;
}

const bclib::DataMatrix& InputAdmixData::getEtaPriorMatrix() const{
  return etaPriorMatrix_;
}

// const bclib::DataMatrix& InputAdmixData::getAlleleFreqMatrix() const {
//     return alleleFreqMatrix_;
// }

// const bclib::DataMatrix& InputAdmixData::getHistoricalAlleleFreqMatrix() const {
//     return historicalAlleleFreqMatrix_;
// }

const bclib::DataMatrix& InputAdmixData::getReportedAncestryMatrix() const{
  return reportedAncestryMatrix_;
}
const Vector_s& InputAdmixData::GetPopLabels() const{
  return HiddenStateLabels;
}

void InputAdmixData::Delete(){
  InputData::Delete();

  //erase string matrices
  for(unsigned i = 0; i < alleleFreqData_.size(); ++i)
    alleleFreqData_[i].clear();
  alleleFreqData_.clear();
   for(unsigned i = 0; i < historicalAlleleFreqData_.size(); ++i)
    historicalAlleleFreqData_[i].clear();
  historicalAlleleFreqData_.clear();
  for(unsigned i = 0; i < etaPriorData_.size(); ++i)
    etaPriorData_[i].clear();
  etaPriorData_.clear();
  for(unsigned i = 0; i < reportedAncestryData_.size(); ++i)
    reportedAncestryData_[i].clear();
  reportedAncestryData_.clear();

  //erase data matrices 
  //alleleFreqMatrix_.clear();
  //historicalAlleleFreqMatrix_.clear();
  etaPriorMatrix_.clear();
  reportedAncestryMatrix_.clear();
}
