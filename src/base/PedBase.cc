//=============================================================================
//
// Copyright (C) 2009  David D. Favro
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
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file PedBase.cc
/// Implementation of the PedBase class.
//=============================================================================

#pragma implementation
#include "PedBase.h"


#include <stdexcept>
#include <typeinfo>
#include <cxxabi.h>


#define ALLOW_UNIMPLEMENTED_X_CHROM	1


static const char NI_PREFIX [] = "Not implemented: ";



//-----------------------------------------------------------------------------
// Helper function: throw nicely-formatted exception.
//-----------------------------------------------------------------------------

// For GNU C, declare the noreturn attribute.
#if defined(__GNUC__)
    static void not_implemented( const char * funcName, const std::type_info & ti )
	__attribute__((noreturn));
#endif

static void not_implemented( const char * funcName, const std::type_info & ti )
    {
    int status;
    std::string message( NI_PREFIX );

    char * const realname = abi::__cxa_demangle( ti.name(), 0, 0, &status );
    message += (status == 0) ? realname : ti.name();
    if ( realname != 0 )
	free( realname );

    message += "::";
    message += funcName;

    throw std::runtime_error( message );
    }



namespace genepi { // ----



PedBase::~PedBase()
    {
    }


void PedBase::DeleteGenotypes() { not_implemented( "DeleteGenotypes()", typeid(*this) ); }

void PedBase::HMMIsBad(bool /*loglikisbad*/) { not_implemented( "HMMIsBad()", typeid(*this) ); }

void PedBase::setOutcome(double*) { not_implemented( "setOutcome()", typeid(*this) ); }

void PedBase::setCovariates(double*) { not_implemented( "setCovariates()", typeid(*this) ); }

void PedBase::setGenotypesToMissing() { not_implemented( "setGenotypesToMissing()", typeid(*this) ); }

void PedBase::SetMissingGenotypes() { not_implemented( "SetMissingGenotypes()", typeid(*this) ); }


const double* PedBase::getAdmixtureProps(bool /* isXChrom */)const
    {
    not_implemented( "getAdmixtureProps()", typeid(*this) );
    }


const std::vector<hapPair > &PedBase::getPossibleHapPairs(unsigned int /*locus*/)const { not_implemented( "getPossibleHapPairs()", typeid(*this) ); }

const int* PedBase::getSampledHapPair(int /*locus*/)const { not_implemented( "getSampledHapPair()", typeid(*this) ); }

bool PedBase::GenotypeIsMissing(unsigned int /*locus*/)const { not_implemented( "GenotypeIsMissing()", typeid(*this) ); }

bool PedBase::simpleGenotypeIsMissing(unsigned /*locus*/)const { not_implemented( "simpleGenotypeIsMissing()", typeid(*this) ); }

bool PedBase::isHaploidatLocus(unsigned /*j*/)const { not_implemented( "isHaploidatLocus()", typeid(*this) ); }

bool PedBase::isHaploidIndividual()const { not_implemented( "isHaploidIndividual()", typeid(*this) ); }

double PedBase::getLogLikelihood(const Options& , const bool /*forceUpdate*/, const bool /*store*/) { not_implemented( "getLogLikelihood()", typeid(*this) ); }

double PedBase::getLogLikelihoodXChr(const Options& , const bool /*forceUpdate*/, const bool /*store*/)
    {
    #if ALLOW_UNIMPLEMENTED_X_CHROM
	return 0.0;
    #else
	not_implemented( "getLogLikelihoodXChr()", typeid(*this) );
    #endif
    }

void PedBase::storeLogLikelihood(const bool /*setHMMAsOK*/) { not_implemented( "storeLogLikelihood()", typeid(*this) ); }

double PedBase::getLogLikelihoodAtPosteriorMeans(const Options& /*options*/) { not_implemented( "getLogLikelihoodAtPosteriorMeans()", typeid(*this) ); }

void PedBase::GetLocusAncestry(int /*locus*/, int /*Ancestry*/[2])const { not_implemented( "GetLocusAncestry()", typeid(*this) ); }

void PedBase::GetLocusAncestry(int /*chrm*/, int /*locus*/, int /*Ancestry*/[2])const { not_implemented( "GetLocusAncestry()", typeid(*this) ); }

int PedBase::GetLocusAncestry(int, int, int)const { not_implemented( "GetLocusAncestry()", typeid(*this) ); }


void PedBase::SampleHiddenStates(const Options& /*options*/) { not_implemented( "SampleHiddenStates()", typeid(*this) ); }


void PedBase::SampleHapPair(unsigned /*chr*/, unsigned /*jj*/, unsigned /*locus*/, AlleleFreqs * /*A*/,
			bool /*skipMissingGenotypes*/, bool /*annealthermo*/, bool /*UpdateCounts*/)
	{ not_implemented( "SampleHapPair()", typeid(*this) ); }


void PedBase::UpdateAlleleCounts(unsigned /*j*/, unsigned /*jj*/, unsigned /*locus*/, AlleleFreqs */*A*/, bool /*annealthermo*/)const { not_implemented( "UpdateAlleleCounts()", typeid(*this) ); }


void PedBase::SampleMissingOutcomes(bclib::DataMatrix * /*Outcome*/, const std::vector<bclib::Regression*>& /*R*/) { not_implemented( "SampleMissingOutcomes()", typeid(*this) ); }


unsigned int PedBase::getMyNumber() const { not_implemented( "getMyNumber()", typeid(*this) ); }

unsigned int PedBase::getIndex() const { not_implemented( "getIndex()", typeid(*this) ); }



// Methods from AdmixedIndividual:

void PedBase::drawInitialAdmixtureProps(const AlphaType & /*alpha*/) { not_implemented( "drawInitialAdmixtureProps()", typeid(*this) ); }


void PedBase::SetGenotypeProbs(int /*j*/, int /*jj*/, unsigned /*locus*/, bool /*chibindicator*/) { not_implemented( "SetGenotypeProbs()", typeid(*this) ); }


void PedBase::AnnealGenotypeProbs(int /*j*/, const double /*coolness*/) { not_implemented( "AnnealGenotypeProbs()", typeid(*this) ); }


void PedBase::ResetSufficientStats() { not_implemented( "ResetSufficientStats()", typeid(*this) ); }


void PedBase::SampleJumpIndicators(bool /*sampleArrivals*/) { not_implemented( "SampleJumpIndicators() ", typeid(*this) ); }


void PedBase::SampleTheta( int /*iteration*/, double * /*SumLogTheta*/, const bclib::DataMatrix * /*Outcome*/,
	const DataType * /*OutcomeType*/, const std::vector<double> & /*lambda*/, int /*NumCovariates*/,
	bclib::DataMatrix * /*Covariates*/, const std::vector<const double*> & /*beta*/, const PopThetaType & /*poptheta*/,
	const AdmixOptions & /*options*/, const AlphaType & /*alpha*/,
	double /*DInvLink*/, double /*dispersion*/, CopyNumberAssocTest & /*ancestryAssocTest*/, bool /*RW*/, bool /*anneal*/) { not_implemented( "SampleTheta()", typeid(*this) ); }


void PedBase::SampleRho(const AdmixOptions& /*options*/, double /*rhoalpha*/,
	double /*rhobeta*/, bool /*updateSumLogRho*/) { not_implemented( "SampleRho()", typeid(*this) ); }


void PedBase::FindPosteriorModes(const AdmixOptions & /*options*/, const AlphaType &/*alpha*/,
	  double /*rhoalpha*/, double /*rhobeta*/, AlleleFreqs* /*A*/, ofstream &/*modefile*/) { not_implemented( "FindPosteriorModes()", typeid(*this) ); }


const RhoType & PedBase::getRho() const { not_implemented( "getRho()", typeid(*this) ); }


void PedBase::resetStepSizeApproximator(int /*k*/) { not_implemented( "resetStepSizeApproximator()", typeid(*this) ); }


void PedBase::UpdateScores(const AdmixOptions& /*options*/, bclib::DataMatrix * /*Outcome*/, bclib::DataMatrix * /*Covariates*/,
    const vector<bclib::Regression*> & /*R*/, AffectedsOnlyTest & /*affectedsOnlyTest*/,
    CopyNumberAssocTest & /*ancestryAssocTest*/ ) { not_implemented( "UpdateScores()", typeid(*this) ); }


void PedBase::setChibNumerator(const AdmixOptions& /*options*/, const AlphaType & /*alpha*/, double /*rhoalpha*/,
	double /*rhobeta*/, chib * /*MargLikelihood*/, AlleleFreqs */*A*/) { not_implemented( "setChibNumerator()", typeid(*this) ); }


void PedBase::updateChib(const AdmixOptions& /*options*/, const AlphaType & /*alpha*/, double /*rhoalpha*/,
	double /*rhobeta*/, chib * /*MargLikelihood*/, AlleleFreqs */*A*/) { not_implemented( "updateChib()", typeid(*this) ); }


double PedBase::getSumrho()const { not_implemented( "getSumrho()", typeid(*this) ); }


double PedBase::getLogLikelihoodOnePop() { not_implemented( "getLogLikelihoodOnePop()", typeid(*this) ); }


double PedBase::getLogPosteriorTheta()const { not_implemented( "getLogPosteriorTheta()", typeid(*this) ); }


double PedBase::getLogPosteriorRho()const { not_implemented( "getLogPosteriorRho()", typeid(*this) ); }


double PedBase::getLogPosteriorAlleleFreqs()const { not_implemented( "getLogPosteriorAlleleFreqs()", typeid(*this) ); }


void PedBase::WritePosteriorMeans(ostream& /*os*/, unsigned /*samples*/, bool /*globalrho*/) const { not_implemented( "WritePosteriorMeans()", typeid(*this) ); }


void PedBase::WritePosteriorMeansXChr(ostream& /*os*/, unsigned /*samples*/) const
    {
    #if ! ALLOW_UNIMPLEMENTED_X_CHROM
	not_implemented( "WritePosteriorMeansXChr()", typeid(*this) );
    #endif
    }


void PedBase::WritePosteriorMeansLoci(ostream& /*os*/) const { not_implemented( "WritePosteriorMeansLoci()", typeid(*this) ); }

void PedBase::setOddsRatios(const genepi::cvector<double>& /*psi*/)
    {
    #if ! ALLOW_UNIMPLEMENTED_X_CHROM
	not_implemented( "setOddsRatios()", typeid(*this) );
    #endif
    }


} // ---- end namespace genepi
