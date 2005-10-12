#include "StratificationTest.h"
#include <numeric>
#include "gsl/gsl_eigen.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

StratificationTest::StratificationTest()
{
   NumberOfTestLoci = 0;
   count = 0;
   T = 0;
}

void StratificationTest::Initialize( AdmixOptions* options, Genome &Loci, Chromosome **Chr, IndividualCollection *IC, LogWriter *Log )
{
  if(options->getStratificationTest() ){
    float DistanceFromLast = 0;
    //OLD CODE
//     for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
//       DistanceFromLast += Loci.GetDistance(j);
//       if( DistanceFromLast > 10 ){ // set cutoff to 10 morgans so that only unlinked loci will be included
// 	// should have a more specific way to code unlinked loci than setting Distance=100
// 	// this algorithm selects the first locus on each chromosome
// 	// preferable to select the most informative locus (greatest heterozygosity)
// 	if( Loci(j)->GetNumberOfStates() == 2 ){ // test uses only diallelic loci
// 	  // should be able to use multi-allelic loci by grouping alleles into 2 bins
// 	  TestLoci.push_back(j);
// 	  NumberOfTestLoci++;
// 	  DistanceFromLast = 0;
// 	}
//       }
//     }
    //NEW CODE
    for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
      int j = -1;
      double max = 0;
      double n1, n2;
      double ExpHet;//expected heterozygosity
      //select most informative locus (with greatest exp heteroxygosity), from those with
      //less than 5% missing genotypes, on each chromosome
      for(unsigned locus = 0; locus < Chr[c]->GetSize(); ++locus){
	if((*Chr[c])(locus)->GetNumberOfStates() == 2){// test uses only diallelic loci
	  //count number of missing genotypes at locus
	  long count = 0;
	  for(int i = 0; i < IC->getSize(); ++i){
	    Individual *ind = IC->getIndividual(i);
	    if(ind->IsMissing(locus)){
	      ++count;
	    }
	  }
	  if( ( (double)count / (double)(IC->getSize()) ) < 0.1){//exclude loci with >10% missing genotypes
	    n1 = (double) GetAlleleCounts(locus, 1, IC);//# copies of allele1
	    n2 = (double) GetAlleleCounts(locus, 2, IC);//# copies of allele2
	    ExpHet = 2.0 * (n1*n2) / ((n1+n2)*(n1+n2)); 
	    if( ExpHet > max){
	      max  = ExpHet;
	      j = Chr[c]->GetLocus(locus);//gets locus number (on genome) of this locus
	    }
	  }
	}
      }
      //j is the most informative locus
      if(j > -1){
	TestLoci.push_back(j);
	NumberOfTestLoci++;
	DistanceFromLast = 0;
      }
    }
    
    
    if( NumberOfTestLoci < 2 ){
      Log->logmsg(true,"Too few unlinked loci to run stratification test\n");
      options->setStratificationTest(false);
    }
    else{
      Log->logmsg(true, NumberOfTestLoci);
      Log->logmsg(true, " loci used in stratification test.\n");
      for(int i = 0; i < NumberOfTestLoci; i++){
	Log->write(Loci(TestLoci[i])->GetLabel(0));Log->write("\n");
      }
      ModelIndicator = options->isRandomMatingModel();
      
      OpenOutputFile(options->getStratTestFilename(),Log);
    }
  }
  else{
    Log->logmsg(true,"No test for residual population stratification.\n");
  }
}

void StratificationTest::calculate( IndividualCollection* individuals, double** AlleleFreqs, vector<vector<int> > ChrmAndLocus, 
				    int Populations )
{
  // matrix of (observed minus expected copies allele 1) scores for each individual at each locus
  gsl_matrix *popX = gsl_matrix_calloc( individuals->getSize(), NumberOfTestLoci ); 
  // matrix of replicate scores
  gsl_matrix *popRepX = gsl_matrix_calloc( individuals->getSize(), NumberOfTestLoci );
  //bool flag = false;
  vector<unsigned short> genotype(2, 0);
  int ancestry[2];
  for( int j = 0; j < NumberOfTestLoci; j++ ){
    int jj = TestLoci[j];
    double* freqs = AlleleFreqs[jj];  // array of length (NumberOfStates-1)*Populations
    for( int i = 0; i < individuals->getSize(); i++ ){
      Individual* ind = individuals->getIndividual(i);
      unsigned short **genotypeArray = ind->getGenotype(jj);
      // recode as vector<unsigned short>
      genotype[0] = genotypeArray[0][0];
      genotype[1] = genotypeArray[0][1];
      ind->GetLocusAncestry( ChrmAndLocus[jj][0], ChrmAndLocus[jj][1], ancestry );
      //if( genotype[0] != genotype[1] ){ // if heterozygous
      // genotype = SampleHeterozygotePhase( freqs, ancestry ); // sample phase conditional on ordered diploid ancestry 
      //} else 
      if( genotype[0] == 0 ){ // if genotype is missing, sample it
	genotype = SimGenotypeConditionalOnAncestry( freqs, ancestry );
      }
      // ProbAllele1 = Prob( allele 1 ) conditional on individual admixture
      vector<double> ProbAllele1 = GenerateExpectedGenotype( ind, freqs, Populations );
      vector<unsigned short> repgenotype = SimGenotypeConditionalOnAdmixture( ProbAllele1 );
      // vector<unsigned short> repgenotype = SimGenotypeConditionalOnAncestry( freqs, ancestry );
      // Calculate score X = Obs0 + Obs1 - Expected0 - Expected1, where Obs0, Obs1 are coded 1 for allele 1, 0 for allele 2
      // Obs0 + Obs1 = 4 - genotype[0] - genotype[1]
      gsl_matrix_set(popX, i, j, (double)(4 - genotype[0] - genotype[1]) - ProbAllele1[0] - ProbAllele1[1]);
      gsl_matrix_set(popRepX, i, j ,(double)(4 - repgenotype[0] - repgenotype[1]) - ProbAllele1[0] - ProbAllele1[1]);
    }  
  }
  
  // covariance matrix for (observed minus expected copies allele 1) scores
  gsl_matrix *Cov = gsl_matrix_calloc(NumberOfTestLoci, NumberOfTestLoci);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, popX, popX, 0, Cov  );//Cov = popX' * popX

  // covariance matrix for replicate scores
  gsl_matrix *RepCov = gsl_matrix_calloc(NumberOfTestLoci, NumberOfTestLoci);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, popRepX, popRepX, 0, RepCov  );//RepCov = popRepX' * popRepX

  //compute eigenvalues
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc (NumberOfTestLoci);//allocates workspace for eigenvalue calculations
  gsl_vector *EigenValues = gsl_vector_calloc(NumberOfTestLoci);
  gsl_eigen_symm (Cov, EigenValues, w);

  gsl_vector *RepEigenValues = gsl_vector_calloc(NumberOfTestLoci);
  gsl_eigen_symm(RepCov, RepEigenValues, w);

  gsl_eigen_symm_free(w); 

  //normalise eigenvalues
  const double Sum = accumulate(EigenValues->data, EigenValues->data + NumberOfTestLoci, 0.0);
  gsl_vector_scale (EigenValues, 1.0/Sum) ;
  const double RepSum = accumulate(RepEigenValues->data, RepEigenValues->data + NumberOfTestLoci, 0.0);
  gsl_vector_scale (RepEigenValues, 1.0/RepSum) ;
  //find max eigenvalues
  double maxEigenValue = gsl_vector_max(EigenValues);
  double maxRepEigenValue = gsl_vector_max(RepEigenValues);
  //output to file
  outputstream << maxEigenValue << "\t" << maxRepEigenValue << endl;
  //increment test statistic
  if( maxEigenValue < maxRepEigenValue ){
    T++;
  }
  count++;
  //clean up
  gsl_vector_free(EigenValues);
  gsl_vector_free(RepEigenValues);
  gsl_matrix_free(Cov);
  gsl_matrix_free(RepCov);
  gsl_matrix_free(popX);
  gsl_matrix_free(popRepX);
}

vector<double>
StratificationTest::GenerateExpectedGenotype( Individual* ind, const double* freqs, const int Populations )
{
  vector<double> pA(2,0);
  for( int k = 0; k < Populations; k++ ){
    pA[0] += freqs[ k ] * ind->getAdmixtureProps()[k];
    if( ModelIndicator )
      pA[1] += freqs[ k ] * ind->getAdmixtureProps()[k + Populations];
    else
      pA[1] += freqs[ k ] * ind->getAdmixtureProps()[k];       
  }
  return pA; // vector of length 2 specifying probs of allele 1 on each gamete
}

vector<unsigned short>
StratificationTest::SimGenotypeConditionalOnAdmixture( const vector<double> ProbAllele1 )
{
  vector<unsigned short> repgenotype(2,0);
  if( ProbAllele1[0] > myrand() )
    repgenotype[0] = 1;
  else
    repgenotype[0] = 2;
  if( ProbAllele1[1] > myrand() )
    repgenotype[1] = 1;
  else
    repgenotype[1] = 2;
  return repgenotype;
}

vector<unsigned short>
StratificationTest::SimGenotypeConditionalOnAncestry( const double* freqs, const int ancestry[2] )
{
  vector<unsigned short> repgenotype(2,0);
  if( freqs[ ancestry[0] ] > myrand() )
    repgenotype[0] = 1;
  else
    repgenotype[0] = 2;
  if( freqs[ ancestry[1] ] > myrand() )
    repgenotype[1] = 1;
  else
    repgenotype[1] = 2;
  return repgenotype;
}

vector<unsigned short> StratificationTest::SampleHeterozygotePhase( const double* freqs, const int Ancestry[2] )
{
  vector<unsigned short> genotype(2, 0);
  double q1 = freqs[ Ancestry[0] ] * ( 1 - freqs[ Ancestry[1] ] );
  double q2 = freqs[ Ancestry[1] ] * ( 1 - freqs[ Ancestry[0] ] );
  if( myrand() > q1 / ( q1 + q2 ) ){
    genotype[0] = 2;
    genotype[1] = 1;
  }
  else{
    genotype[0] = 1;
    genotype[1] = 2;
  }
  return genotype;
}

void StratificationTest::OpenOutputFile( const char * OutputFilename, LogWriter *Log){

  outputstream.open(OutputFilename, ios::out );
  if( !outputstream ){
    Log->logmsg(false,"ERROR: Couldn't open stratificationtestfile");
    Log->logmsg(false,OutputFilename);
    Log->logmsg(false,"\n");
    exit( 1 );
  }
  Log->logmsg(true,"Writing results of test for residual population stratification to ");
  Log->logmsg(true, OutputFilename);
  Log->logmsg(true,"\n");
  outputstream << "T_obs" << "\t" << "T_rep\n";
}

void StratificationTest::Output(LogWriter *Log){
  Log->logmsg(true, "\nStratification test: posterior predictive check probability "); 
  Log->logmsg(true, (float)T/count);
  Log->logmsg(true,"\n\n");
}

//this function should really be in CompositeLocus or AlleleFreqs or IndividualCollection
int StratificationTest::GetAlleleCounts(int locus, int a, IndividualCollection *IC)
{
  /**
   * returns a count of the copies of allele a at a comp locus
   * only works for diallelic loci. 
   */

  int AlleleCounts = 0;

  for(int i = 0; i < IC->getSize(); ++i){
    Individual *ind = IC->getIndividual(i);
    if(!ind->IsMissing(locus)){
      unsigned short **genotype = ind->getGenotype(locus);
	if(genotype[0][0] == a){
	  AlleleCounts++;
	}
	if(genotype[0][1] == a){
	  AlleleCounts++;
	}
    }
  }
  return AlleleCounts;
}
