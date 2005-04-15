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
  StationaryDist = new double[D];
  //Likelihood.SetNumberOfElementsWithDimensions( L, D, 1);
  GenotypeProbs = new double *[L];
  for(int i=0;i<L;++i)GenotypeProbs[i] = new double[D];
  CodedStates = new int[L];
  Tpat = new double*[populations] ;
  Tmat = new double*[populations];
  for(int i=0;i<populations;++i){
    Tpat[i] = new double[populations];
    Tmat[i] = new double[populations];
  }
  _product1 = new double[populations];
  _product2 = new double[populations];

}

void
Chromosome::ResetStuffForX()
{

  D = populations;
  SampleStates.SetDimensions( L, D );
  if(StationaryDist != NULL) delete StationaryDist;
  StationaryDist = new double[ D ];
  //Likelihood.SetNumberOfElementsWithDimensions( L, D, 1);
  for(int i=0;i<L;++i){
    delete GenotypeProbs[i];
    GenotypeProbs[i] = new double[D];
  }

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
  delete CodedStates;
  delete StationaryDist;
  delete[] GenotypeProbs;
  delete[] Tpat;
  delete[] Tmat;

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
Chromosome::UpdateParameters(Individual* ind, AlleleFreqs *A, Matrix_d& Admixture, AdmixOptions* options, Vector_d f[],
			     bool fixedallelefreqs, bool diploid )
//Obtains stationary distribution and transition probs for HMM and updates forward and backward probabilities
//Admixture - matrix of admixture proportions
//f - 
{

  int locus, d;

  if(diploid){
    int Mcol = options->getModelIndicator(); //which col to use as maternal ancestry, col 1 for a randommatingmodel
    //col 0 is paternal ancestry
    // Construct stationary distribution
    d = 0;
    for( int k = 0; k < populations; k++ ){
      for( int kk = 0; kk < populations; kk++ ){
	StationaryDist[ d ] = Admixture(k,0) * Admixture(kk,Mcol);
	d++;
      }
    }
    
    if(L > 1){
      locus = GetLocus( 0 );
      for( int j = 0; j < L - 1; j++ ){
	// Construct Haploid transition matrices, Tpat and Tmat, then use them to construct transition matrix
	
	for(int i=0; i<populations; i++){
	  _product1[i] = f[0](locus+1) * Admixture(i,0);
	  _product2[i] = f[1](locus+1) * Admixture(i,Mcol);
	}
	
	for(int i=0; i<populations; i++){
	  for(int j=0; j<i; j++){
	    Tpat[i][j] = Admixture(j,0) - _product1[j];
	    Tmat[i][j] = Admixture(j,Mcol) - _product2[j];
	    Tpat[j][i] = Admixture(i,0) - _product1[i];
	    Tmat[j][i] = Admixture(i,Mcol) - _product2[i];
	  }
	  Tpat[i][i] = Admixture(i,0) + f[0](locus+1) - _product1[i];
	  Tmat[i][i] = Admixture(i,Mcol) + f[1](locus+1) - _product2[i];
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
		//TransitionProbs(j)( row, col ) = Tpat[k1][kk1]*Tmat[k2][kk2];
		SampleStates.SetTProb(j, row, col, Tpat[k1][kk1]*Tmat[k2][kk2] );
		col++;
	      }
	    }
	    row++;
	  }
	}
	
	locus++;
      }
    }
    //this could be neater
    //case of a chromosome with a single locus - set transition probs to zero
    else for(int j= 0; j< L-1; j++)for(int k1 = 0; k1<populations*populations;k1++)for(int k2 = 0; k2<populations*populations; k2++)
      SampleStates.SetTProb(j, k1, k2, 0.0);
  }

  else{//haploid
    Vector_d mu = Admixture.GetColumn(0);
    for( int j = 0; j < D; j++ )
      StationaryDist[j] = mu(j); 
    
    if(L >1){
      locus = GetLocus( 0 );
      //construct Transition Probs  
      for( int j = 0; j < L - 1; j++ ){
	int H = mu.GetNumberOfElements();
	//TransitionProbs(j).SetNumberOfElements(H, H);
	
	double _product[H];
	for(int i=0; i<H; i++){
	  _product[i] = f[0](locus+1) * mu(i);
	}
	
	for(int i=0; i<H; i++){
	  for(int jj=0; jj<i; jj++){
	    //TransitionProbs(j)(i,jj) = mu(jj) - _product[jj];
	    SampleStates.SetTProb(j,i,jj, mu(jj) - _product[jj]);
	  }
	  //TransitionProbs(j)(i,i) = mu(i) + f[0](locus+1) - _product[i];
	  SampleStates.SetTProb(j,i,i, mu(i) + f[0](locus+1) - _product[i]);
	  for(int jj=i+1; jj<H; jj++){
	    //TransitionProbs(j)(i,jj) = mu(jj) - _product[jj];
	    SampleStates.SetTProb(j,i,jj, mu(jj) - _product[jj]);
	  }
	}
	
	locus++;
      }
    }
    //this could be neater
    //case of a chromosome with a single locus - set transition probs to zero
    else for(int j= 0; j< L-1; j++)for(int k1 = 0; k1<populations*populations;k1++)for(int k2 = 0; k2<populations*populations; k2++)
      SampleStates.SetTProb(j, k1, k2, 0.0);
  }
 
  // Construct array of probs of genotype given ancestry states
  locus = GetLocus( 0 );
  for( int j = 0; j < L; j++ ){
    d = 0;
    if( !(ind->IsMissing(locus)) ){
       A->GetGenotypeProbs(&Prob, locus, ind->getGenotype(locus), ind->getPossibleHaplotypes(locus), true, fixedallelefreqs );
      for( int k = 0; k < populations; k++ ){
	for( int kk = 0; kk < populations; kk++ ){
	  GenotypeProbs[j][d] = Prob( k, kk );
	  d++;
	}
      }
    }
    else{
      for( int k = 0; k < D; k++ ) GenotypeProbs[j][k] = 1.0;
    }
    locus++;
  }

  bool test = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());

  SampleStates.UpdateFwrdBckwdProbabilities( StationaryDist, GenotypeProbs, test );

}

void Chromosome::SampleForLocusAncestry(Matrix_i *OrderedStates, bool isdiploid)
//Ordered States - pointer to ((2 or 1) x size) int matrix to store the locus ancestry states
//isdiploid - indicator for diploid (true) or haploid (false)
{

  // Sample    
  SampleStates.Sample(CodedStates);
  for( int j = 0; j < L; j++ ){
    if(isdiploid){
      //     OrderedStates( 0, j ) = 0;
      (*OrderedStates)( 0, j ) = (int)(CodedStates[j] / populations);
      (*OrderedStates)( 1, j ) = (CodedStates[j] % populations);
    }
    else //haploid
      {
	(*OrderedStates)(0, j ) = CodedStates[j];
      }
  }
   
}

void Chromosome::getAncestryProbs( int j, Matrix_d *AncestryProbs ){
  //sets conditional probabilities of ancestry at locus j
  //One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
  //i.e. (i,2) = p_{ii}
  //     (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}
  //     (i,0) = 1.0 - (i,1) - (i,2)
  //where p's are probs in StateProbs
  
  AncestryProbs->SetNumberOfElements(populations, 3);
  double *StateProbs;
  StateProbs = new double[D];//possibly should keep this at class scope
  
   SampleStates.GetStateProbs(StateProbs, j);
  
  for( int k1 = 0; k1 < populations; k1++ ){
    (*AncestryProbs)(k1,2) = StateProbs[ ( populations + 1 ) * k1 ];
    for( int k2 = 0 ; k2 < populations; k2++ )
      (*AncestryProbs)(k1,1) += StateProbs[k1*populations +k2] + StateProbs[k2*populations +k1];
    (*AncestryProbs)(k1,1)-= 2.0*(*AncestryProbs)(k1,2);
    (*AncestryProbs)(k1,0) = 1.0 - (*AncestryProbs)(k1,1) - (*AncestryProbs)(k1,2);
  }
  delete StateProbs;
}


double Chromosome::getLogLikelihood()
{
   return SampleStates.getLikelihood();
}
