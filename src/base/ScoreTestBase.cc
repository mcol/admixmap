#include "ScoreTestBase.h"
#include <sstream>
#include <numeric>
#include "gsl/gsl_cdf.h"
#include <math.h>
#include "bcppcl/linalg.h"

using namespace std;

ScoreTestBase::~ScoreTestBase(){
  if(outputfile.is_open())outputfile.close();
}

void ScoreTestBase::OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname, bool Robj){
  outputstream->open(filename, ios::out);
  if(!outputstream->is_open()){
    string error_string = "ERROR: could not open ";
    error_string.append(filename);
    throw(error_string);
  }
  Log << testname << " written to " << filename << "\n";
  if(Robj)
    //start writing R object
    *outputstream << "structure(.Data=c(";

  onFirstLine = true;
}


///generic scalar score test
void ScoreTestBase::OutputScalarScoreTest( int iterations, ofstream* outputstream, string label,
					const double score, const double scoresq, const double info, bool final)
{
  string sep = final? "\t" : ",";
  const double Score = score / ( double )iterations;
  const double CompleteInfo = info / (double) iterations;
  const double MissingInfo = scoresq / (double) iterations - Score * Score;
  const double ObservedInfo = CompleteInfo - MissingInfo;

  if(!onFirstLine){
    //if not the first line, output a separator after previous output
    *outputstream << sep;
  }
  //now start a new line (for ease of human reading)
  *outputstream << endl;
  //after this function call there will be at least one line in output
  onFirstLine = false;

  //output label
  *outputstream << "\"" << label << "\"" << sep;
  if(final)
    *outputstream << double2R(Score, 3)        << sep
		  << double2R(CompleteInfo, 3) << sep
		  << double2R(ObservedInfo, 3) << sep;

  if(CompleteInfo > 0.0 && (MissingInfo < CompleteInfo)) {
    const double PercentInfo = 100.0 * ObservedInfo / CompleteInfo;
    if(final) *outputstream << double2R(PercentInfo, 2) << sep;

    if(!final || PercentInfo > 10.0){ //only output p-values in final table if >10% info extracted
      const double zscore = Score / sqrt( ObservedInfo );
      const double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
      if(final){
	*outputstream << double2R(zscore,3) << sep 
		      << double2R(pvalue);
      }
      else //not final table - output log p-value
	*outputstream << double2R(-log10(pvalue));
    }
    else //final table and %info is <10
      *outputstream << "NA" << sep << "NA";

  }// negative CI or MissingInfo > CompleteInfo
  else{
    if(final)*outputstream << "NA" << sep << "NA";
    *outputstream << "NA" ;
  }
  //*outputstream << endl;
}

void ScoreTestBase::OutputRaoBlackwellizedScoreTest( ofstream* outputstream, string label,
						     const double score, const double scoresq, const double varscore, 
						     const double info, bool final )
{
  if(numUpdates == 0){
    throw string("Unable to output scoretest as no updates have been made");
  }

  string separator = final? "\t" : ",";

  if(!onFirstLine){
    //if not the first line, output a separator after previous output
    *outputstream << separator;
  }
  //now start a new line (for ease of human reading)
  *outputstream << endl;
  //after this function call there will be at least one line in output
  onFirstLine = false;


  *outputstream << "\"" << label << "\"" << separator;

  const double scaleFactor = 1.0 / (double)numUpdates;      
  const double EU = score * scaleFactor;
  const double VU = varscore * scaleFactor;
  const double missing = scoresq * scaleFactor - EU * EU + VU;
  const double complete =  info * scaleFactor;
  
  if(final){
    *outputstream << double2R(EU, 3)                                << separator//score
		  << double2R(complete, 3)                          << separator//complete info
		  << double2R(complete - missing, 3)                << separator;//observed info

  }
  if(complete > 0.0){
    const double PercentInfo = 100.0*(complete - missing)/complete;
    if(final){
      *outputstream << double2R(PercentInfo, 2) << separator//%observed info
		    << double2R(100.0*(VU/complete), 2)                 << separator//%missing info attributable to locus ancestry
		    << double2R(100.0*(missing-VU)/complete, 2)         << separator;//%remainder of missing info      
    }
    if(missing < complete) {
      if(!final || PercentInfo > 10.0){ //only output p-values in final table if >10% info extracted
      const double zscore = EU / sqrt( complete - missing );
      const double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
      if(final){
	  *outputstream << double2R(zscore,3) << separator 
			<< double2R(pvalue) ;
	}
      else //not final table - output log p-value
	*outputstream << double2R(-log10(pvalue));
      }
      else //final table and %info is <10
	*outputstream << "NA" << separator << "NA";
    }
    else{// MissingInfo > CompleteInfo
      if(final)*outputstream << "NA" << separator;
      *outputstream << "NA";
    }
  }
  else{//complete info <= 0
    if(final)
      *outputstream << "NA" << separator << "NA" << separator << "NA" << separator << "NA" << separator;
    *outputstream << "NA" ; 
  }
  //*outputstream << endl;
}

///generic vector score test
void ScoreTestBase::OutputScoreTest( int iterations, ofstream* outputstream, unsigned dim, vector<string> labels,
				  const double* score, const double* scoresq, const double* info, bool final, unsigned dim2)
{
  //given cumulative scores, square of scores and info, of dimension dim, over iterations, computes expectation of score, complete info and observed info and outputs to output stream along with a summary chi-square statistic and p-value. Also performs scalar test for each element.
  //if final=false, only the log(-pvalue)'s are printed

  double *ScoreVector = 0, *CompleteInfo = 0, *ObservedInfo = 0;
  string sep = final? "\t" : ",";

  if(!onFirstLine){
    //if not the first line, output a separator after previous output
    *outputstream << sep;
  }
  //now start a new line (for ease of human reading)
  *outputstream << endl;
  //after this function call there will be at least one line in output
  onFirstLine = false;

  
  ScoreVector = new double[dim];
  copy(score, score+dim, ScoreVector);
  scale_matrix(ScoreVector, 1.0/( iterations), dim, 1);
  
  CompleteInfo = new double[dim*dim];
  copy(info, info + dim*dim, CompleteInfo);
  scale_matrix(CompleteInfo, 1.0/( iterations), dim, dim);
  
  ObservedInfo = new double[dim*dim];
  for(unsigned d1 = 0; d1 < dim; ++d1)for(unsigned d2 = 0; d2 < dim; ++d2)
    ObservedInfo[d1*dim + d2] = CompleteInfo[d1*dim+d2] + ScoreVector[d1]*ScoreVector[d2] -
      scoresq[d1*dim+d2]/( iterations );
  for( unsigned k = 0; k < dim; k++ ){
    // ** output labels
    *outputstream << labels[k] << sep;
    if(final){
      *outputstream  << double2R(ScoreVector[k], 3) << sep
		     << double2R(CompleteInfo[k*dim+k], 3) << sep//prints diagonal of CI matrix
		     << double2R(ObservedInfo[k*dim+k], 3) << sep;//   "      "     "  MI   "
      if(CompleteInfo[k*dim+k]>0.0)
	*outputstream<< double2R(100*ObservedInfo[k*dim+k] / CompleteInfo[k*dim+k], 2) << sep;//%Observed Info
      else 
	*outputstream << "NA" << sep;
    }
    double zscore = ScoreVector[ k ] / sqrt( ObservedInfo[k*dim+k] );
    if(final)*outputstream  << double2R(zscore, 3) << sep;//z-score
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final)*outputstream << double2R(pvalue) << sep;
    else *outputstream << double2R(-log10(pvalue));
    // if not last allele at locus, output unquoted "NA" in chi-square column
    if( final && k != dim - 1 ){
      *outputstream  << "NA";
    }
  }//end loop over alleles
  if(final){
    double chisq=0.0;
    try{
      if(dim2==dim) chisq = GaussianQuadraticForm(ScoreVector, ObservedInfo, dim);
      else chisq = GaussianMarginalQuadraticForm( dim2, ScoreVector, ObservedInfo, dim );//marginalise over first dim2 elements
      if(chisq < 0.0)
	*outputstream << "NA";
      else *outputstream << double2R(chisq);
    }
    catch(...){//in case ObservedInfo is rank deficient
      *outputstream  << "NA";
    }
  }
  //TODO:?? output p-value for chisq
	
  delete[] ScoreVector;
  delete[] CompleteInfo;
  delete[] ObservedInfo;
}

///finishes writing scoretest output as 3-dimensional R object
void ScoreTestBase::R_output3DarrayDimensions(ofstream* stream, const vector<int> dim, const vector<string> labels)
{
  *stream << endl << ")," << endl;
  *stream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    *stream << dim[i];
    if(i != dim.size() - 1){
      *stream << ",";
    }
  }
  *stream << ")," << endl;

  //output dimnames for the first dimension (column labels)
  *stream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    *stream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      *stream << ",";
    }
  }
  //output 'chcracter(0)' for last 2 dimensions
  *stream << "), character(0), character(0)))" << endl;
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
