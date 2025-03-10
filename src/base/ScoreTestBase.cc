//=============================================================================
//
// Copyright (C) 2006, 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file ScoreTestBase.cc
/// Implementation of the ScoreTestBase class.
//=============================================================================

#include "ScoreTestBase.h"
#include "bclib/linalg.h"
#include "bclib/DelimitedFileWriter.h"
#include "gsl/gsl_cdf.h"
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace std;

ScoreTestBase::ScoreTestBase(){
  test = false;
}
ScoreTestBase::~ScoreTestBase(){
}

/// Generic scalar score test
void ScoreTestBase::OutputScalarScoreTest( int iterations, bclib::DelimitedFileWriter& outputstream, string label,
					const double score, const double scoresq, const double info, bool final)
{
  const double Score = score / ( double )iterations;
  const double CompleteInfo = info / (double) iterations;
  const double MissingInfo = scoresq / (double) iterations - Score * Score;
  const double ObservedInfo = CompleteInfo - MissingInfo;

  //output label
  outputstream << ("\"" + label + "\"") ;
  if(final)
    outputstream << double2R(Score, 3)  
		  << double2R(CompleteInfo, 3)
		  << double2R(ObservedInfo, 3);

  if(CompleteInfo > 0.0 && ObservedInfo > 0.0) {
    const double PercentInfo = 100.0 * ObservedInfo / CompleteInfo;
    if(final) outputstream << double2R(PercentInfo, 2);

    if(PercentInfo >= 10.0){ //only output p-values if >10% info extracted
      const double zscore = Score / sqrt( ObservedInfo );
      const double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
      if(final){
	//output zscore and pvalue in final table
	outputstream << double2R(zscore,3) 
		      << double2R(pvalue) ;
      }
      else //not final table - output log p-value
	outputstream << double2R(-log10(pvalue));
    }
    else{ // %info is <10
      if(final)//2 NAs in final table (z-score, p-value)
	outputstream << "NA" << "NA";
      else//1 NA in cumulative output (log p-value)
	outputstream << "NA";
    }
  }

  // negative CompleteInfo or MissingInfo > CompleteInfo
  else {
    if (final)
      outputstream << "NA" << "NA"; // %Info,  z-score
    //NA for (log)p-value
    outputstream << "NA" ;
  }
  outputstream << bclib::newline;
}

void ScoreTestBase::OutputRaoBlackwellizedScoreTest(  bclib::DelimitedFileWriter& outputstream, string label,
						     const double score, const double scoresq, const double varscore, 
						     const double info, bool final )
{
  if(numUpdates == 0){
    throw string("Unable to output scoretest as no updates have been made");
  }

  outputstream << ("\"" + label + "\"");

  const double scaleFactor = 1.0 / (double)numUpdates;      
  const double EU = score * scaleFactor;
  const double VU = varscore * scaleFactor;
  const double missing = scoresq * scaleFactor - EU * EU + VU;
  const double complete =  info * scaleFactor;
  const double observed = complete - missing;

  if (final)
    outputstream << double2R(EU, 3)        // score
                 << double2R(complete, 3)  // complete info
                 << double2R(observed, 3); // observed info

  if(complete > 0.0){

    // observed info (%)
    const double PercentInfo = 100.0 * observed / complete;

    // missing info attributable to hidden state (%)
    const double Missing1 = 100.0 * VU / complete;

    // remainder of missing info (%)
    const double Missing2 = 100.0 * (missing - VU) /complete;

    if (final) {
      outputstream << (PercentInfo >= 0.0 ? double2R(PercentInfo, 2) : "NA");
      outputstream << double2R(Missing1, 2) << double2R(Missing2, 2);
    }
    if(missing < complete) {
      if(PercentInfo >= 10.0){ //only output p-values if >10% info extracted
	const double zscore = EU / sqrt(observed);
	const double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
	if (final) //output zscore and pvalue in final table
	  outputstream << double2R(zscore, 3) << double2R(pvalue);
	else //not final table - output log p-value
	  outputstream << double2R(-log10(pvalue));
      }
      else{ // %info is <10
	if(final)//2 NAs in final table
	  outputstream << "NA" << "NA";
	else//1 NA in cumulative output
	  outputstream << "NA";
      }
    }
    else{// MissingInfo > CompleteInfo
      if(final)outputstream << "NA";
      outputstream << "NA";
    }
  }

  // complete info <= 0
  else {
    if(final)
      outputstream << "NA" << "NA" << "NA" << "NA" ;
    outputstream << "NA" ; 
  }
  outputstream << bclib::newline;
}

/// Generic vector score test
void ScoreTestBase::OutputScoreTest( int iterations,  bclib::DelimitedFileWriter& outputstream, unsigned dim, vector<string> labels,
				  const double* score, const double* scoresq, const double* info, bool final, unsigned dim2)
{
  //given cumulative scores, square of scores and info, of dimension dim, over iterations, computes expectation of score, complete info and observed info and outputs to output stream along with a summary chi-square statistic and p-value. Also performs scalar test for each element.
  //if final=false, only the log(-pvalue)'s are printed

  double *ScoreVector = 0, *CompleteInfo = 0, *ObservedInfo = 0;

  ScoreVector = new double[dim];
  copy(score, score+dim, ScoreVector);
  bclib::scale_matrix(ScoreVector, 1.0/( iterations), dim, 1);
  
  CompleteInfo = new double[dim*dim];
  copy(info, info + dim*dim, CompleteInfo);
  bclib::scale_matrix(CompleteInfo, 1.0/( iterations), dim, dim);
  
  ObservedInfo = new double[dim*dim];
  for(unsigned d1 = 0; d1 < dim; ++d1)for(unsigned d2 = 0; d2 < dim; ++d2)
    ObservedInfo[d1*dim + d2] = CompleteInfo[d1*dim+d2] + ScoreVector[d1]*ScoreVector[d2] -
      scoresq[d1*dim+d2]/( iterations );
  for( unsigned k = 0; k < dim; k++ ){
    // ** output labels
    outputstream << labels[k] ;
    if(final){
      outputstream  << double2R(ScoreVector[k], 3) 
		     << double2R(CompleteInfo[k*dim+k], 3) //prints diagonal of CI matrix
		     << double2R(ObservedInfo[k*dim+k], 3);//   "      "     "  MI   "
      if(CompleteInfo[k*dim+k]>0.0)
	outputstream<< double2R(100*ObservedInfo[k*dim+k] / CompleteInfo[k*dim+k], 2);//%Observed Info
      else 
	outputstream << "NA" ;
    }
    double zscore = ScoreVector[ k ] / sqrt( ObservedInfo[k*dim+k] );
    if(final)outputstream  << double2R(zscore, 3) ;//z-score
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final){
      outputstream << double2R(pvalue) ;
      if( k != dim - 1 ){
	outputstream  << "NA" << bclib::newline;
      }
     }
    else outputstream << double2R(-log10(pvalue)) << bclib::newline;
    // if not last allele at locus, output unquoted "NA" in chi-square column
  }//end loop over alleles
  if(final){
    double chisq=0.0;
    try{
      if(dim2==dim) chisq = bclib::GaussianQuadraticForm(ScoreVector, ObservedInfo, dim);
      else chisq = bclib::GaussianMarginalQuadraticForm( dim2, ScoreVector, ObservedInfo, dim );//marginalise over first dim2 elements
      if(chisq < 0.0)
	outputstream << "NA";
      else outputstream << double2R(chisq);
    }
    catch(...){//in case ObservedInfo is rank deficient
      outputstream  << "NA";
    }
    outputstream << bclib::newline;
  }
  //TODO:?? output p-value for chisq
	
  delete[] ScoreVector;
  delete[] CompleteInfo;
  delete[] ObservedInfo;
}

///converts a double to a string for R to read
//useful only for converting infs and nans to "NaN"
string ScoreTestBase::double2R( double x )
{
  if( isnan(x) || isinf(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
string ScoreTestBase::double2R( double x, int precision )
{
  if( isnan(x) || isinf(x) )
    return "NaN";
  else{
    stringstream ret;

    if(x < numeric_limits<float>::max( ))//in case x too large to represent as floating point number
      ret << setiosflags(ios::fixed) << setprecision(precision);
    ret << x;
    return( ret.str() );
  }
}
