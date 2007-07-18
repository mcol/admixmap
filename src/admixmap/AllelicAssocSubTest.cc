#include "AllelicAssocSubTest.h"
#include "CompositeLocus.h"
#include "ScoreTestBase.h"
#include "bclib/linalg.h"
#include "bclib/FileWriter.h"
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <string>
#include <sstream>

using namespace::std;

unsigned AllelicAssocSubTest::NumCovars;

AllelicAssocSubTest::AllelicAssocSubTest(){
  Score = 0;
  Info = 0;
  dim = 0;
}

AllelicAssocSubTest::AllelicAssocSubTest(unsigned d) :
  dim(d){

  Score = new double[ dim + NumCovars ];
  Info = new double[( dim + NumCovars) * (dim + NumCovars )];
}

void AllelicAssocSubTest::SetNumCovars(unsigned n){
  NumCovars = n;
}

AllelicAssocSubTest::~AllelicAssocSubTest(){
  delete[] Score;
  delete[] Info;
}

unsigned AllelicAssocSubTest::getDim()const{
  return dim;
}

void AllelicAssocSubTest::Reset(){
  fill(Score, Score+dim+NumCovars, 0.0);
  fill(Info, Info+(dim+NumCovars)*(dim+NumCovars), 0.0);
}

void AllelicAssocSubTest::Update(const vector<int> Counts, const double* covariates,  
				 double YMinusEY, double phi, double DInvLink)
{
  double* x = new double[ NumCovars + dim ];

  // ** Set x co-ordinate for regression parameter under test
  copy(Counts.begin(), Counts.end(), x);
 
  // ** Set x-co-ordinates of covariates in model
  x[ dim ] = 1.0;//intercept
  for( unsigned k = 1; k < NumCovars; k++ )
    x[ dim+k ] = covariates[k];
  //if(!options->getHapMixModelIndicator())
  //copy(covariates+1, covariates+NumCovars, x+dim+1);

  //accumulate score and info
  for( unsigned k = 0; k < NumCovars + dim; k++ ){
    Score[ k ] += phi * x[k] * YMinusEY;
    for( unsigned kk = 0; kk < NumCovars + dim; kk++ )
      Info[ k*(NumCovars+dim) + kk ] += x[ k ] * x[ kk ] * phi*DInvLink;
  }
  delete[] x;
}

// corrects for covariance between score and regression parameters
// then accumulates sumscore, sumscoresq and suminfo
void AllelicAssocSubTest::CentreAndSum(unsigned dim, double *score, double* info, 
				       double *sumscore, double* sumscoresq, double* suminfo)
{
  double *cscore = new double[dim];//centred score
  double *cinfo = new double[dim*dim];//centred info

  bclib::CentredGaussianConditional( dim, score, info, cscore, cinfo, NumCovars+dim );
  for(unsigned d = 0; d < dim; ++d){
    *(sumscore + d) += cscore[d];
    for(unsigned dd = 0; dd < dim; ++dd){
      *(sumscoresq + d*dim +dd) += cscore[d] * cscore[dd];
      *(suminfo + d*dim +dd) += cinfo[d*dim+dd];
    }
  }
  delete[] cscore;
  delete[] cinfo;
}

// next few function calculate score tests from the cumulative sums of
// the score, score squared, and information score and info can be
// scalars, or respectively a vector and a matrix

// generic scalar score test
//TODO: move output of NA in chisq column outside as it is only required if along with vector tests
void AllelicAssocSubTest::OutputScalarScoreTest(bclib::FileWriter& outputstream, string label, 
						const double score, const double scoresq, 
						const double info, bool final, unsigned numUpdates)
{
  double Score = score / (double) numUpdates;
  double CompleteInfo = info / (double) numUpdates;
  double MissingInfo = scoresq / ( double)numUpdates - Score * Score;
  double ObservedInfo = CompleteInfo - MissingInfo;

  //output label
  const string quotedlabel = "\"" + label + "\"" ;
  outputstream << quotedlabel;

  if(final)
    outputstream << ScoreTestBase::double2R(Score, 3)       
		  << ScoreTestBase::double2R(CompleteInfo, 3)
		  << ScoreTestBase::double2R(ObservedInfo, 3);

  if( (MissingInfo < CompleteInfo) && (CompleteInfo > 0.0) ) {
    double PercentInfo = 100.0*ObservedInfo / CompleteInfo;
    double zscore = Score / sqrt( ObservedInfo );
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final)
      outputstream << ScoreTestBase::double2R(PercentInfo, 2)
		    << ScoreTestBase::double2R(zscore,3)
		    << ScoreTestBase::double2R(pvalue);
    else
      outputstream << ScoreTestBase::double2R(-log10(pvalue));
  }
  else{
    if(final)outputstream << "NA" << "NA" ;
    outputstream << "NA" ;
  }
  if(final)outputstream << "NA";//NA in chisquare column in final table 
  outputstream << bclib::newline;
}

//generic vector score test
void AllelicAssocSubTest::OutputScoreTest( bclib::FileWriter& outputstream, unsigned dim, vector<string> labels,
					   const double* score, const double* scoresq, const double* info, 
					   bool final, unsigned dim2, unsigned numUpdates)
{
  //given cumulative scores, square of scores and info, of dimension dim, over iterations, computes expectation of score, complete info and observed info and outputs to output stream along with a summary chi-square statistic and p-value. Also performs scalar test for each element.
  //if final=false, only the log(-pvalue)'s are printed

  double *ScoreVector = 0, *CompleteInfo = 0, *ObservedInfo = 0;

  ScoreVector = new double[dim];
  copy(score, score+dim, ScoreVector);
  bclib::scale_matrix(ScoreVector, 1.0/( numUpdates), dim, 1);
  
  CompleteInfo = new double[dim*dim];
  copy(info, info + dim*dim, CompleteInfo);
  bclib::scale_matrix(CompleteInfo, 1.0/(double) numUpdates, dim, dim);
  
  ObservedInfo = new double[dim*dim];
  for(unsigned d1 = 0; d1 < dim; ++d1)for(unsigned d2 = 0; d2 < dim; ++d2)
    ObservedInfo[d1*dim + d2] = CompleteInfo[d1*dim+d2] + ScoreVector[d1]*ScoreVector[d2] -
      scoresq[d1*dim+d2]/(double) numUpdates;
  for( unsigned k = 0; k < dim; k++ ){

    // ** output labels
    outputstream << labels[k];
    if(final){
      outputstream  << ScoreTestBase::double2R(ScoreVector[k], 3)
		     << ScoreTestBase::double2R(CompleteInfo[k*dim+k], 3) //prints diagonal of CI matrix
		     << ScoreTestBase::double2R(ObservedInfo[k*dim+k], 3);//   "      "     "  MI   "
      if(CompleteInfo[k*dim+k]>0.0)
	outputstream<< ScoreTestBase::double2R(100.0*ObservedInfo[k*dim+k] / CompleteInfo[k*dim+k], 2) ;//%Observed Info
      else 
	outputstream << "NA" ;
    }
    double zscore = ScoreVector[ k ] / sqrt( ObservedInfo[k*dim+k] );
    if(final)outputstream  << ScoreTestBase::double2R(zscore, 3);//z-score
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final){
      outputstream << ScoreTestBase::double2R(pvalue);
      if( k != dim - 1 ){
	outputstream  << "NA" << bclib::newline;
      }
    }
    else outputstream << ScoreTestBase::double2R(-log10(pvalue)) << bclib::newline;
    // if not last allele at locus, output unquoted "NA" in chi-square column
  }//end loop over alleles
  if(final){
    double chisq=0.0;
    try{
      if(dim2==dim) chisq = bclib:: GaussianQuadraticForm(ScoreVector, ObservedInfo, dim);
      else chisq = bclib::GaussianMarginalQuadraticForm( dim2, ScoreVector, ObservedInfo, dim );//marginalise over first dim2 elements
      if(chisq < 0.0)
	outputstream << "NA" ;
      else outputstream << ScoreTestBase::double2R(chisq) ;
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

/*************************************
           SNP Test
*************************************/

SNPTest::SNPTest() : AllelicAssocSubTest(1){
  SumScore = 0.0;
  SumScore2 = 0.0;
  SumInfo = 0.0;
}

SNPTest::~SNPTest(){
}

void SNPTest::Update(const int* const happair, CompositeLocus* const, const double* covariates,  
		    double YMinusEY, double phi, double DInvLink){
  vector<int> allele2counts;
  allele2counts.push_back( (happair[0]==1) + (happair[1]==1) );
  // Locus->getAlleleCounts(2, happair)[0];
  
  AllelicAssocSubTest::Update(allele2counts, covariates, YMinusEY, phi, DInvLink);
}

void SNPTest::Accumulate(){
  CentreAndSum(dim, Score, Info, &SumScore, &SumScore2, &SumInfo); 
}

void SNPTest::Output(bclib::FileWriter& outfile, const CompositeLocus* const Locus, 
		     bool final,unsigned numUpdates){
  OutputScalarScoreTest(outfile, Locus->GetLabel(0), SumScore, SumScore2, 
			SumInfo, final, numUpdates);
}

/*************************************
           Multiallelic locus test
*************************************/
MultiAllelicLocusTest::MultiAllelicLocusTest(unsigned d) : AllelicAssocSubTest(d){
  //next two lines may not be necessary as these arrays are sized later

  SumScore = new double[ dim ];
  SumInfo = new double[dim * dim];
  SumScore2 = new double[dim * dim];
  fill(SumScore, SumScore+dim, 0.0);
  fill(SumInfo, SumInfo + dim*dim, 0.0);
  fill(SumScore2, SumScore2 + dim*dim, 0.0);
}

MultiAllelicLocusTest::~MultiAllelicLocusTest(){
  delete[] SumScore;
  delete[] SumScore2;
  delete[] SumInfo;
}

void MultiAllelicLocusTest::Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
				   double YMinusEY, double phi, double DInvLink){
  // count alleles / haplotypes      
  vector<int> counts;
  for( unsigned k = 0; k < dim; k++ ){
    counts.push_back(Locus->getAlleleCounts(k+1, happair)[0] );
  }

  AllelicAssocSubTest::Update(counts, covariates, YMinusEY, phi, DInvLink);
}

void MultiAllelicLocusTest::Accumulate(){
  CentreAndSum(dim, Score, Info, SumScore, SumScore2, SumInfo); 
}

void MultiAllelicLocusTest::Output(bclib::FileWriter& outfile, const CompositeLocus* const Locus, 
				   bool final, unsigned numUpdates){

  vector<string> labels;
  for(unsigned i = 0; i < dim; ++i){
    stringstream ss;
    ss << "\"" << Locus->GetLabel(0) << "(" << i+1 << ")\"";
    labels.push_back(ss.str());
  }
  OutputScoreTest(outfile, dim, labels, 
		  SumScore, SumScore2, SumInfo, final, dim, numUpdates);
}

/*************************************
           Haplotype Test
*************************************/
HaplotypeTest::HaplotypeTest(unsigned d) 
  : MultiAllelicLocusTest(d){
}

void HaplotypeTest::Resize(unsigned d){
  dim = d;

  delete[] Score;
  delete[] Info;
  Score = new double[ dim + NumCovars ];
  Info = new double[( dim + NumCovars) * (dim + NumCovars )];

  delete[] SumScore;
  delete[] SumScore2;
  delete[] SumInfo;
  SumScore = new double[ dim ];
  SumInfo = new double[dim * dim];
  SumScore2 = new double[dim * dim];
  fill(SumScore, SumScore+dim, 0.0);
  fill(SumInfo, SumInfo + dim*dim, 0.0);
  fill(SumScore2, SumScore2 + dim*dim, 0.0);
}

HaplotypeTest::~HaplotypeTest(){
}

void HaplotypeTest::Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
				double YMinusEY, double phi, double DInvLink){
  //count numbers of each haplotype
  vector<int> counts = Locus->getHaplotypeCounts(happair);
  
  AllelicAssocSubTest::Update(counts, covariates, YMinusEY, phi, DInvLink);
}

void HaplotypeTest::Accumulate(){
    CentreAndSum(dim, Score, Info, SumScore, SumScore2, SumInfo); 
}

void HaplotypeTest::Output(bclib::FileWriter& outfile, const CompositeLocus* const Locus, 
			   bool final, unsigned numUpdates){
  //here, dim = NumberOfMergedHaplotypes
  //create labels as "locuslabel","haplabel"
  vector<string> labels;
  const int *hap = 0;
  const string sep = final ? "\t" : ",";//separator
  const unsigned NumberOfLoci = Locus->GetNumberOfLoci();

  for(unsigned i = 0; i < dim; ++i){
    stringstream ss;
    ss  << "\"" << Locus->GetLabel(0) << "\""<< sep;
    if( i < dim - 1 ){
      hap = Locus->GetHapLabels(i);
      ss  << "\"";
      for( unsigned kk = 0; kk < NumberOfLoci - 1; kk++ ){
	ss  << hap[kk] << "-";
      }
      ss  << hap[NumberOfLoci - 1] << "\"";
    }
    else
      ss  << "\"others\"";
    labels.push_back(ss.str());
  }
  OutputScoreTest(outfile, dim, labels, SumScore, SumScore2, SumInfo, final, dim-1, numUpdates);
}

/*************************************
           Within-haplotype Test
*************************************/
WithinHaplotypeTest::WithinHaplotypeTest(unsigned d){
  dim = d;

  WithinHaplotypeScore = new double*[dim];
  WithinHaplotypeInfo = new double*[dim];
  for(unsigned jj = 0; jj < dim; ++jj){
    WithinHaplotypeScore[jj] = new double[1 + NumCovars];
    WithinHaplotypeInfo[jj] = new double[(1 + NumCovars)*(1 + NumCovars)];
  }

  SumWithinHaplotypeScore  = new double[ dim ];
  SumWithinHaplotypeScore2 = new double[ dim ];
  SumWithinHaplotypeInfo   = new double[ dim ];
  fill(SumWithinHaplotypeScore, SumWithinHaplotypeScore+dim, 0.0);
  fill(SumWithinHaplotypeScore2, SumWithinHaplotypeScore2+dim, 0.0);
  fill(SumWithinHaplotypeInfo, SumWithinHaplotypeInfo+dim, 0.0);
}

WithinHaplotypeTest::~WithinHaplotypeTest(){
  for(unsigned jj = 0; jj < dim; ++jj){
    delete[] WithinHaplotypeScore[jj];
    delete[] WithinHaplotypeInfo[jj];
  }

  delete[] WithinHaplotypeScore;
  delete[] WithinHaplotypeInfo;
  delete[] SumWithinHaplotypeScore;
  delete[] SumWithinHaplotypeScore2;
  delete[] WithinHaplotypeInfo;
}

void WithinHaplotypeTest::Reset(){
  for(unsigned jj = 0; jj < dim; ++jj){
    fill(WithinHaplotypeScore[jj], WithinHaplotypeScore[jj]+NumCovars+1, 0.0);
    fill(WithinHaplotypeInfo[jj], WithinHaplotypeInfo[jj]+(NumCovars+1)*(NumCovars+1), 0.0);
  }
}

void WithinHaplotypeTest::Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
				double YMinusEY, double phi, double DInvLink){
  //update score and info for each simple locus within a compound locus
  // 	for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
  // 	  vector<int> a(1, allele2Counts[0]);
  // 	  UpdateAlleleScores(ScoreWithinHaplotype[locus][l], InfoWithinHaplotype[locus][l], ind->getAdmixtureProps(), a, 
  // 			     YMinusEY, phi, DInvLink);
  // 	}

  //count copies of allele2
  const vector<int> allele2Counts = Locus->getAlleleCounts(2, happair);
  UpdateWithinHaplotypeAssociationTest(covariates, allele2Counts, YMinusEY,phi , DInvLink);

}

// This function calculates score for allelic association at each simple locus within a compound locus
void WithinHaplotypeTest::UpdateWithinHaplotypeAssociationTest( const double* covariates, 
							  const vector<int> allele2Counts, 
							  double YMinusEY, double phi, double DInvLink)
{
  double* x = new double[ NumCovars + 1 ];

  x[ NumCovars ] = 1.0;
  for( unsigned k = 0; k < NumCovars - 1; k++ )
    x[ k + 1 ] = covariates[k];//?? should be [k+1]

  for( unsigned l = 0; l < dim; l++ ){//loop over simple loci within compound locus
    x[0] = (double)allele2Counts[l];

    for( unsigned k = 0; k < NumCovars + 1; k++ ){
      WithinHaplotypeScore[l][ k ] += phi * x[k] * YMinusEY;
      for( unsigned kk = 0; kk < NumCovars + 1; kk++ )
	WithinHaplotypeInfo[l][ k*(NumCovars+1) + kk ] += phi*DInvLink * x[ k ] * x[ kk ];
    }
  }
  delete[] x;
}

void WithinHaplotypeTest::Accumulate(){
  for( unsigned l = 0; l < dim; l++ ){
    CentreAndSum(1, WithinHaplotypeScore[l], WithinHaplotypeInfo[l], &(SumWithinHaplotypeScore[ l ]),
		 &(SumWithinHaplotypeScore2[l]), &(SumWithinHaplotypeInfo[ l ]));
  }
}

void WithinHaplotypeTest::Output(bclib::FileWriter& outfile, const CompositeLocus* const Locus, 
			   bool final, unsigned numUpdates){

  //for(int i = 0; i < (*Lociptr)(j)->GetNumberOfLoci(); ++i)labels.push_back("\""+(*Lociptr)(j)->GetLabel(i)+"\"");
  for(unsigned simplelocus = 0; simplelocus < dim; ++simplelocus){
    string simplelocuslabel = Locus->GetLabel(simplelocus);
    OutputScalarScoreTest(outfile, simplelocuslabel, SumWithinHaplotypeScore[simplelocus], 
			  SumWithinHaplotypeScore2[simplelocus], SumWithinHaplotypeInfo[simplelocus], 
			  final, numUpdates);
  }
}
