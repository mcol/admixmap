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
      if( DistanceFromLast > 0.05 ){
	if( Loci(j)->GetNumberOfStates() == 2 ){
	  TestLoci.push_back(j);
	  NumberOfTestLoci++;
	  DistanceFromLast = 0;
	}
      }
    }
    if( NumberOfTestLoci < 2 ){
       Log->logmsg(true,"Can't run stratification test with this data set.\n");
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

void StratificationTest::calculate( IndividualCollection* individuals, Matrix_d *AlleleFreqs, vector<vector<int> > ChrmAndLocus )
{
  Matrix_d popX( individuals->getSize(), NumberOfTestLoci );
  Matrix_d popRepX( individuals->getSize(), NumberOfTestLoci );
  bool flag = false;

  for( int j = 0; j < NumberOfTestLoci; j++ ){
    int jj = TestLoci[j];
    Matrix_d freqs = AlleleFreqs[jj];
    for( int i = 0; i < individuals->getSize(); i++ ){
      Individual* ind = individuals->getIndividual(i);
     unsigned short **genotype = ind->getGenotype(jj);
      if( genotype[0][0] ){
	Vector_i ancestry = ind->GetLocusAncestry( ChrmAndLocus[jj][0], ChrmAndLocus[jj][1] );
	vector<unsigned short >repgenotype = GenerateRepGenotype( freqs, ancestry );
	vector<double> pA = GenerateExpectedGenotype( ind, freqs );
	if( genotype[0][0] != genotype[0][1] ){
	  genotype = SampleForOrderedSNiP( freqs, ancestry );
	  flag = true;//need to delete genotype
	}
// pA = P( allele 1 ).
// Calculate error X = Y - pA, where Y = 1 for allele 1 and Y = 2 for allele 2.
// The following code actually calculates X / Xrep = -Y + pA.
        popRepX( i, j ) = (double)(repgenotype[0]+repgenotype[1])-4+pA[0]+pA[1];
        popX( i, j ) = (double)(genotype[0][0]+genotype[0][1])-4+pA[0]+pA[1];
	//if(flag){delete[] genotype[0];delete[] genotype;}
      }
    }
  }
  
  Vector_d EigenValues;
  Vector_d RepEigenValues;
  Matrix_d Cov;
  Matrix_d RepCov;

  Cov = popX.Transpose() * popX;
  RepCov = popRepX.Transpose() * popRepX;
  Cov.Eigenvalue2( &EigenValues );
  RepCov.Eigenvalue2( &RepEigenValues );
  EigenValues /= EigenValues.Sum();
  RepEigenValues /= RepEigenValues.Sum();
  if( EigenValues.MaximumElement() < RepEigenValues.MaximumElement() ){
     T++;
  }
  count++;
}

vector<double>
StratificationTest::GenerateExpectedGenotype( Individual* ind, const Matrix_d& freqs )
{
  vector<double> pA(2,0);
  for( int k = 0; k < freqs.GetNumberOfCols(); k++ ){
    pA[0] += freqs( 0, k ) * ind->getAdmixtureProps()[k];
    if( ModelIndicator )
      pA[1] += freqs( 0, k ) * ind->getAdmixtureProps()[k + freqs.GetNumberOfCols()];
    else
       pA[1] += freqs( 0, k ) * ind->getAdmixtureProps()[k];       
  }
  return pA;
}

vector<unsigned short>
StratificationTest::GenerateRepGenotype( const Matrix_d& freqs, const Vector_i& ancestry )
{
  vector<unsigned short> repgenotype(2,0);
  if( freqs( 0, ancestry(0) ) > myrand() )
    repgenotype[0] = 1;
  else
    repgenotype[0] = 2;
  if( freqs( 0, ancestry(1) ) > myrand() )
    repgenotype[1] = 1;
  else
    repgenotype[1] = 2;
  return repgenotype;
}

unsigned short **StratificationTest::SampleForOrderedSNiP( const Matrix_d& freqs, const Vector_i& Ancestry )
{
  unsigned short **genotype;
  genotype = new unsigned short*[1];
  genotype[0] = new unsigned short[2];
  
  double q1 = freqs( 0, Ancestry(0) ) * ( 1 - freqs( 0, Ancestry(1) ) );
  double q2 = freqs( 0, Ancestry(1) ) * ( 1 - freqs( 0, Ancestry(0) ) );
  if( myrand() > q1 / ( q1 + q2 ) ){
    genotype[0][0] = 2;
    genotype[0][1] = 1;
  }
  else{
    genotype[0][0] = 1;
    genotype[0][1] = 2;
  }
  return genotype;
}

// float StratificationTest::getStatistic()
// {
//    return (float)T/count;
// }

void StratificationTest::OpenOutputFile( const char * OutputFilename, LogWriter *Log){

  outputstream.open(OutputFilename, ios::out );
  if( !outputstream ){
    Log->logmsg(false,"ERROR: Couldn't open stratificationtestfile");
    Log->logmsg(false,OutputFilename);
    Log->logmsg(false,"\n");
    exit( 1 );
  }
  Log->logmsg(true,"Writing results of test for residual population stratification to ");
  Log->logmsg(true,OutputFilename);
  Log->logmsg(true,"\n");
  //write header
  outputstream << "Post. Pred. Check Prob." <<endl;
}

void StratificationTest::Output(){
   outputstream << (float)T/count << endl;
}
