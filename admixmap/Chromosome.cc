#include "Chromosome.h"
#include "Individual.h"
#include <iostream>

#define PR(x) cerr << #x << " = " << x << endl;

using namespace std;

Chromosome::Chromosome(int size, int start, int inpopulations) : Genome(size)
{
  _startLoci = start;
  populations = inpopulations;
  D = populations * populations;
  L = GetNumberOfCompositeLoci();
  SampleStates.SetDimensions( L, D );
  StationaryDist.SetNumberOfElements( D, 1 );
  // as an array of transition matrices is part of the definition of an HMM, 
  // it would make more sense to have one HMM object for each chromosome, 
  // and for this HMM object to store transition probs
  TransitionProbs.SetNumberOfElementsWithDimensions( L - 1, D, D );
  Likelihood.SetNumberOfElementsWithDimensions( L, D, 1 );
  //Tpat and Tmat are the paternal and maternal haploid transition probability matrices 
  //Tpat = 0;
  //Tmat = 0;
}

void
Chromosome::ResetStuffForX()
{
  D = populations;
  SampleStates.SetDimensions( L, D );
  StationaryDist.SetNumberOfElements( D, 1 );
  TransitionProbs.SetNumberOfElementsWithDimensions( L - 1, D, D );
  Likelihood.SetNumberOfElementsWithDimensions( L, D, 1 );
}

void
Chromosome::SetLabel( int, string label )
{
   _Label = label;
}

string
Chromosome::GetLabel( int )
{
   return _Label;
}

Chromosome::~Chromosome()
{

}

//
// Returns the number of the num'th compositelocus on this chromosome
// eg if chromosome 2 consists of compositeloci 5,6,7 and 8,
// GetLocus(i) returns 5 + i
int
Chromosome::GetLocus(int num){
  return _startLoci + num;
}

int
Chromosome::GetSize(){
  return GetNumberOfCompositeLoci();
}

void
Chromosome::UpdateParameters(Individual* ind, AlleleFreqs *A, Matrix_d& _ancestry, AdmixOptions* options, vector< Vector_d >& f,
			     bool fixedallelefreqs )
{
 
  Matrix_d Prob;
  MatrixArray_d empty;
  // Construct stationary distribution
  Vector_d mu_2, mu_1 = _ancestry.GetColumn(0); //maternal
  if( options->getModelIndicator() )
    mu_2 = _ancestry.GetColumn(1); //paternal
  else
    mu_2 = _ancestry.GetColumn(0); //paternal
  
  int d = 0;
  for( int k = 0; k < populations; k++ ){
    for( int kk = 0; kk < populations; kk++ ){
      StationaryDist( d, 0 ) = mu_1(k) * mu_2(kk);
      d++;
    }
  }
  //  StationaryDist = mu_2.ColumnMatrix();
  
  int locus = GetLocus( 0 );
  for( int j = 0; j < L - 1; j++ ){
    // Construct Haploid transition matrices, Tpat and Tmat, then use them to construct transition matrix

    double **Tpat, **Tmat;
    unsigned int H1, H2;
    H1 = mu_1.GetNumberOfElements();
    H2 = mu_2.GetNumberOfElements();  
    Tpat = new double*[H1];
    Tmat = new double*[H2]; 
 
    double _product1[H1], _product2[H2];
    for(unsigned int i=0; i<H1; i++){
      Tpat[i] = new double[H1];
      _product1[i] = f[0](locus+1) * mu_1(i);
    }

    for(unsigned int i=0; i<H1; i++){
      for(unsigned int j=0; j<i; j++){
	Tpat[i][j] = mu_1(j) - _product1[j];
      }
      Tpat[i][i] = mu_1(i) + f[0](locus+1) - _product1[i];
      for(unsigned int j=i+1; j<H1; j++){
	Tpat[i][j] = mu_1(j) - _product1[j];
      }
    }

    
    for(unsigned int i=0; i<H2; i++){
      Tmat[i] = new double[H2];
      _product2[i] = f[1](locus+1) * mu_2(i);
    }
    
    for(unsigned int i=0; i<H2; i++){
      for(unsigned int j=0; j<i; j++){
	Tmat[i][j] = mu_2(j) - _product2[j];
      }
      Tmat[i][i] = mu_2(i) + f[1](locus+1) - _product2[i];
      for(unsigned int j=i+1; j<H2; j++){
	Tmat[i][j] = mu_2(j) - _product2[j];
      }
    }


    // set elements of diploid transition matrix
    int row = 0;
    for( int k1 = 0; k1 < populations; k1++ ){
      for( int k2 = 0; k2 < populations; k2++ ){
	//           int row = k2 + k1 * populations;
	int col = 0;
	for( int kk1 = 0; kk1 < populations; kk1++ ){
	  for( int kk2 = 0; kk2 < populations; kk2++ ){
	    //                 int col = kk2 + kk1 * populations;
	    TransitionProbs(j)( row, col ) = Tpat[k1][kk1]*Tmat[k2][kk2];
	    col++;
	  }
	}
	row++;
      }
    }
   for(unsigned int i=0; i<H1; i++){
      delete Tpat[i];
      }
   for(unsigned int i=0; i<H2; i++){
      delete Tmat[i];
      }

    delete Tpat;
    delete Tmat;

    locus++;
  }
 
  // Construct likelihood (vector of probs of genotype given ancestry states)
  locus = GetLocus( 0 );
  for( int j = 0; j < L; j++ ){
    d = 0;
    if( ind->IsMissing(locus)[0] != 0 ){
      //vector<unsigned int> genotype = ind->getGenotype(locus);
      //Prob = (*this)(j)->GetLikelihood( genotype, true, fixedallelefreqs );
      Prob = A->GetLikelihood(locus, ind->getGenotype(locus), ind->getPossibleHaplotypes(locus), true, fixedallelefreqs );
      for( int k = 0; k < populations; k++ ){
	for( int kk = 0; kk < populations; kk++ ){
	  Likelihood(j)( d, 0 ) = Prob( k, kk );
	  d++;
	}
      }
    }
    else
      Likelihood(j).SetElements( 1.0 );
    locus++;
  }

  bool test = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  if( L > 1 )
    SampleStates.UpdateParameters( StationaryDist, TransitionProbs,
				   Likelihood, test );
  else
    SampleStates.UpdateParameters( StationaryDist, empty,
				   Likelihood, test );

}

void
Chromosome::UpdateParametersHaploid(Individual* ind, AlleleFreqs *A, Matrix_d& _ancestry, AdmixOptions* options, vector< Vector_d >& f,
bool fixedallelefreqs )
{
  Matrix_d Prob;
  MatrixArray_d empty;
  Vector_d mu = _ancestry.GetColumn(0);
  StationaryDist = mu.ColumnMatrix();
  
  int locus = GetLocus( 0 );
  for( int j = 0; j < L - 1; j++ ){
    //HaploidTransitionMatrix Tmat( mu, f[0](locus+1) );

     int H = mu.GetNumberOfElements();
     TransitionProbs(j).SetNumberOfElements(H, H);
     
     double _product[H];
    for(int i=0; i<H; i++){
      _product[i] = f[0](locus+1) * mu(i);
    }

    for(int i=0; i<H; i++){
      for(int jj=0; jj<i; jj++){
	TransitionProbs(j)(i,jj) = mu(jj) - _product[jj];
      }
      TransitionProbs(j)(i,i) = mu(i) + f[0](locus+1) - _product[i];
      for(int jj=i+1; jj<H; jj++){
	TransitionProbs(j)(i,jj) = mu(jj) - _product[jj];
      }
    }

    //TransitionProbs(j) = Tmat.toMatrix();
     locus++;
  }
  // Construct likelihood
  locus = GetLocus( 0 );
  for( int j = 0; j < L; j++ ){
     if( ind->IsMissing(locus)[0] != 0 ){
       //vector<unsigned int> genotype = ind->getGenotype(locus);
        //Likelihood(j) = (*this)(j)->GetLikelihood( genotype, false, fixedallelefreqs );
       Likelihood(j) = A->GetLikelihood(locus, ind->getGenotype(locus), ind->getPossibleHaplotypes(locus), false, fixedallelefreqs);
     }
     else
        Likelihood(j).SetElements( 1.0 );
     locus++;
  }

  bool test = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  if( L > 1 )
     SampleStates.UpdateParameters( StationaryDist, TransitionProbs,
                                    Likelihood, test );
  else
     SampleStates.UpdateParameters( StationaryDist, empty,
                                    Likelihood, test );
}

Matrix_i
Chromosome::SampleForLocusAncestry(Individual* ind, AlleleFreqs *A)
{
  // D - number of diploid ancestry states
  Vector_i CodedStates;
  Matrix_i OrderedStates(2,L);
  
  // Sample    
  CodedStates = SampleStates.Sample(TransitionProbs);
  for( int j = 0; j < L; j++ ){
    //     OrderedStates( 0, j ) = 0;
    OrderedStates( 0, j ) = (int)(CodedStates(j) / populations);
    OrderedStates( 1, j ) = (CodedStates(j) % populations);
  }
  
  // Update stats for allele freqs
  for( int j = 0; j < L; j++ ){
    int locus = GetLocus( j );
    if( ind->IsMissing(locus)[0] != 0 ){
      //vector<unsigned int> genotype = ind->getGenotype(locus);
      //(*this)(j)->UpdateAlleleCounts( genotype, OrderedStates.GetColumn(j) );
      // why should function UpdateAlleleCounts need to get genotypes if it already has the possible haplotypes 
      // compatible with the genotype? 
      A->UpdateAlleleCounts( locus, ind->getGenotype(locus), ind->getPossibleHaplotypes(locus), OrderedStates.GetColumn(j) );
    }
  }
  
  return( OrderedStates );
}

Vector_i
Chromosome::SampleForHaploidLocusAncestry(Individual* ind, AlleleFreqs* A)
{
  Vector_i OrderedStates;
  OrderedStates = SampleStates.Sample(TransitionProbs);

  for( int j = 0; j < L; j++ ){
     int locus = GetLocus( j );
     if( ind->IsMissing(locus)[0] != 0 ){
       //vector<unsigned int> genotype = ind->getGenotype(locus);
	//(*this)(j)->UpdateAlleleCounts_HaploidData( genotype, OrderedStates(j) );
	A->UpdateAlleleCounts_HaploidData( locus, ind->getGenotype(locus), OrderedStates(j) );
     }
  }

  return( OrderedStates );
}

Matrix_d Chromosome::getExpectedAncestry( int j )
{
  
  //One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
  //i.e. (i,2) = p_{ii}
  //     (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}
  //     (i,0) = 1.0 - (i,1) - (i,2)
  //where p's are probs in StateProbs
  
  Matrix_d AncestryProbs( populations, 3 );
  Vector_d StateProbs;
  
  StateProbs = SampleStates.GetStateProbs(j);//vector of ordered ancestry state probs
  
  for( int k1 = 0; k1 < populations; k1++ ){
    AncestryProbs(k1,2) = StateProbs( ( populations + 1 ) * k1 );
    for( int k2 = 0 ; k2 < populations; k2++ )
      AncestryProbs(k1,1) += StateProbs(k1*populations +k2) + StateProbs(k2*populations +k1);
    AncestryProbs(k1,1)-= 2.0*AncestryProbs(k1,2);
    AncestryProbs(k1,0) = 1.0 - AncestryProbs(k1,1) - AncestryProbs(k1,2);
  }
  
  return AncestryProbs;
}

double Chromosome::getLogLikelihood()
{
   return SampleStates.getLikelihood();
}
