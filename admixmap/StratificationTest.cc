#include "StratificationTest.h"

using namespace std;

StratificationTest::StratificationTest()
{
   NumberOfTestLoci = 0;
   count = 0;
   T = 0;
}

void StratificationTest::Initialize( AdmixOptions* options, Genome &Loci, LogWriter *Log )
{
  if(options->getStratificationTest() ){
    float DistanceFromLast = 0;
    for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
      DistanceFromLast += Loci.GetDistance(j);
      if( DistanceFromLast > 10 ){ // set cutoff to 10 morgans so that only unlinked loci will be included
	// should have a more specific way to code unlinked loci than setting Distance=100
	// this algorithm selects the first locus on each chromosome
	// preferable to select the most informative locus (greatest heterozygosity)
	if( Loci(j)->GetNumberOfStates() == 2 ){ // test uses only diallelic loci
	  // should be able to use multi-allelic loci by grouping alleles into 2 bins
	  TestLoci.push_back(j);
	  NumberOfTestLoci++;
	  DistanceFromLast = 0;
	}
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
       
       OpenOutputFile(options->getDICoutputFilename(),Log);
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
  Matrix_d popX( individuals->getSize(), NumberOfTestLoci ); 
  // matrix of replicate scores
  Matrix_d popRepX( individuals->getSize(), NumberOfTestLoci );
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
      if( genotype[0] != genotype[1] ){ // if heterozygous
	genotype = SampleHeterozygotePhase( freqs, ancestry ); // sample phase conditional on ordered diploid ancestry 
      } else 
	if( genotype[0] == 0 ){ // if genotype is missing, sample it
	  genotype = SimGenotypeConditionalOnAncestry( freqs, ancestry );
	}
      // ProbAllele1 = Prob( allele 1 ) conditional on individual admixture
      vector<double> ProbAllele1 = GenerateExpectedGenotype( ind, freqs, Populations );
      vector<unsigned short> repgenotype = SimGenotypeConditionalOnAdmixture( ProbAllele1 );
      // vector<unsigned short> repgenotype = SimGenotypeConditionalOnAncestry( freqs, ancestry );
      // Calculate score X = Obs0 + Obs1 - Expected0 - Expected1, where Obs0, Obs1 are coded 1 for allele 1, 0 for allele 2
      // Obs0 + Obs1 = 4 - genotype[0] - genotype[1]
      popX( i, j )    = (double)(4 - genotype[0] - genotype[1]) - ProbAllele1[0] - ProbAllele1[1];
      popRepX( i, j ) = (double)(4 - repgenotype[0] - repgenotype[1]) - ProbAllele1[0] - ProbAllele1[1];
    }  
  }
  
  Vector_d EigenValues;
  Vector_d RepEigenValues;
  Matrix_d Cov;  // covariance matrix for (observed minus expected copies allele 1) scores
  Matrix_d RepCov; // covariance matrix for replicate scores
  
  Cov = popX.Transpose() * popX;
  RepCov = popRepX.Transpose() * popRepX;
  Cov.Eigenvalue2( &EigenValues );
  RepCov.Eigenvalue2( &RepEigenValues );
  EigenValues /= EigenValues.Sum();
  RepEigenValues /= RepEigenValues.Sum();
  outputstream << EigenValues.MaximumElement() << "\t" << RepEigenValues.MaximumElement() << endl;
  if( EigenValues.MaximumElement() < RepEigenValues.MaximumElement() ){
    T++;
  }
  count++;
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
  Log->logmsg(true, "Posterior predictive check probability "); 
  Log->logmsg(true, (float)T/count);
  Log->logmsg(true,"\n");
}
