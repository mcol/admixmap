#include "CompositeLocus.h"


using namespace std;

/**
 * NAME
 *
 *   CompositeLocus - represents a composite locus.
 *
 *
 * SYNOPSIS
 *
 *   #include "CompositeLocus.h"
 *
 *   CompositeLocus loc123;
 *   loc123.SetNumberOfLoci(3);
 *
 *   [!! NEED MORE EXAMPLES !!]
 */

/**
 * No-argument constructor for a CompositeLocus object, which represents a
 * composite locus. By default, constructs a simple diallelic locus with fixed 
 * allele frequencies.
 */
CompositeLocus::CompositeLocus()
{
  NumberOfLoci = 1;
  NumberOfAlleles.SetNumberOfElements( NumberOfLoci );
  NumberOfAlleles(0) = 2;
  NumberOfStates = 2;

  
  Vector_i      null_Vector_i(1);
  Vector_d      null_Vector_d(1);
  Matrix_i      null_Matrix_i(1,1);
  Matrix_d      null_Matrix_d(1,1);
  MatrixArray_d null_MatrixArray_d;
  VectorLoop    null_VectorLoop(null_Vector_i);
  
  // not previously defined
  Populations = 0;
  NumberOfMergedHaplotypes = 0;
  NumberOfAlleles = null_Vector_i;
  HaplotypeProbs = null_MatrixArray_d;

  DipLoop = null_VectorLoop;
  HapLoop = null_VectorLoop;
  // std::string Label;
  ScoreGene = null_Matrix_d;
  InfoGene = null_Matrix_d;
  SumScoreGene = null_Matrix_d;
  SumScoreGeneSq = null_Matrix_d;
  SumInfoGene = null_Matrix_d;
  pmult = null_Vector_i;
  MergeHaplotypes = null_Vector_i;

  SumNewScore = null_MatrixArray_d;
  SumNewInfo = null_MatrixArray_d;
  SumNewScoreSq = null_MatrixArray_d;
  HapLabels = null_Matrix_i;
}

CompositeLocus::~CompositeLocus()
{
   delete [] Label;
}

/**
 * Changes number of loci in this composite locus, then sets each
 * locus to be diallelic (as for a SNP).
 *
 * NewNumberOfLoci - the number of loci to be represented by this
 * object.
 */
void CompositeLocus::SetNumberOfLoci( int NewNumberOfLoci )
{
   NumberOfLoci = NewNumberOfLoci;
   NumberOfAlleles.SetNumberOfElements( NumberOfLoci );
   NumberOfAlleles.SetElements( 2 );
   NumberOfStates = (int)pow( 2.0, NumberOfLoci );
}

/**
 * Gets the number of loci in this composite locus.
 * 
 * Returns:
 * the number of loci
 */
int CompositeLocus::GetNumberOfLoci()
{
   return( NumberOfLoci );
}

/**
 * Returns the number of states which can exist within this composite 
 * locus. If this composite locus is composed of a single locus, the 
 * number of states will be equal to the number of alleles. If there 
 * are more than one locus in this composite, the number of states 
 * will be equal to the number of possible haplotypes.
 *
 * Returns:
 * the number of states
 */
int CompositeLocus::GetNumberOfStates()
{
   return( NumberOfStates );
}

void CompositeLocus::SetNumberOfStates( int newNumberOfStates )
{
   NumberOfStates= newNumberOfStates;
}

void CompositeLocus::SetNumberOfPopulations(int pops){
  Populations = pops;
}
/**
 * Sets the number of alleles at a given locus.
 * Exits with an error if the locus doesn't exist.
 *
 * alleles - the number of alleles that exist at a given locus
 * locus - the given locus
 */
void CompositeLocus::SetNumberOfAllelesOfLocus( int locus, int alleles )
{
   if( locus > NumberOfLoci - 1 ){
      cout << "Input to SetNumberOfAllelesOfLocus > NumberOfLoci\n";
      exit(0);
   }

   NumberOfStates /= NumberOfAlleles( locus );
   NumberOfAlleles( locus ) = alleles;
   NumberOfStates *= alleles;
}

/**
 * Gets the number of alleles at a given locus.
 * Exits with an error if the locus does not exist.
 *
 * locus - the locus
 * 
 * returns:
 * the number of alleles at the given locus
 */
int CompositeLocus::GetNumberOfAllelesOfLocus( int locus )
{
   if( locus > NumberOfLoci - 1 ){
      cout << "Input to GetNumberOfAllelesOfLocus > NumberOfLoci\n";
      exit(0);}

   return( NumberOfAlleles( locus ) );
}

/**
 * Extends the composite locus by adding one locus (containing a given
 * number of alleles) to the end of the composite locus.
 *
 * alleles - the number of alleles in the locus to be added
 */

void CompositeLocus::AddLocus( int alleles )
{ 
   NumberOfAlleles.AddElement( NumberOfLoci );
   NumberOfAlleles( NumberOfLoci ) = alleles;
   NumberOfLoci++;
   NumberOfStates *= alleles;
}

void CompositeLocus::SetNumberOfLabels()
{
   Label = new string[ NumberOfLoci ];
}

void CompositeLocus::InitialiseScoreTest(int Populations )
{

   ScoreGene.SetNumberOfElements( Populations, 1 );
   SumScoreGene.SetNumberOfElements( Populations, 1 );
   InfoGene.SetNumberOfElements( Populations, Populations );
   SumInfoGene.SetNumberOfElements( Populations, Populations );
   SumScoreGeneSq.SetNumberOfElements( Populations, Populations );
   SumNewScore.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, 1 );
   SumNewInfo.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, NumberOfStates - 1 );
   SumNewScoreSq.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, NumberOfStates - 1 );
}

void CompositeLocus::InitialiseHaplotypes(Matrix_d &Freqs){
    ConstructHaplotypeProbs(Freqs);
    HaplotypeProbsMAP = HaplotypeProbs;
    SetNoMergeHaplotypes();
 }

/**
 * Sets the name of this composite locus (usually from the 
 * allelefreqs.txt file
 *
 * newlabel - the name of this composite locus
 */
void CompositeLocus::SetLabel( int index, string newlabel )
{
   Label[index] = newlabel;
}

/**
 * Gets the name of this composite locus
 *
 * returns:
 * the name of this composite locus
 */
string CompositeLocus::GetLabel(int index)
{
   return( Label[index] );
}

/**
 * Given the genotype and ancestry of mother and father, returns a
 * haplotype pair.
 * 
 * genotype - a two-element vector of paternal and maternal genotypes
 *   in decimal notation (e.g. 1121 , 1221 ).
 * 
 * Haplotypes - a vector of possible haplotypes compatible with genotype
 *
 * ancestry - a two-element vector of parental ancestry (e.g. 1,0 
 *   might represent european paternal and african maternal).
 *
 * AlleleFreqs - a Matrix of allele frequencies, with one column per population
 *
 * returns:
 * two-element vector containing the haplotype pair
 */
Vector_i CompositeLocus::SampleHaplotypePair(Vector_i Haplotypes, Vector_i ancestry)
{
   Vector_i hap(2);
   //   if( NumberOfLoci > 1 ){
      int i;
      Vector_d Probs;
   
      Probs.SetNumberOfElements( 2*Haplotypes.GetNumberOfElements() );
      for( int k = 0; k < Haplotypes.GetNumberOfElements(); k++ ){
         Probs( 2*k )     = HaplotypeProbs( Haplotypes(k) )( ancestry(0), ancestry(1) );
         Probs( 2*k + 1 ) = HaplotypeProbs( Haplotypes(k) )( ancestry(1), ancestry(0) );
      }
      i = SampleFromDiscrete2( &Probs );
      DipLoop.Reset();
      DipLoop.Increment( Haplotypes( i/2 ) );
      if( i % 2 ){
         hap(0) = (DipLoop.GetCountParent( 0 )*pmult).Sum();
         hap(1) = (DipLoop.GetCountParent( 1 )*pmult).Sum();
      }
      else{
         hap(1) = (DipLoop.GetCountParent( 0 )*pmult).Sum();
         hap(0) = (DipLoop.GetCountParent( 1 )*pmult).Sum();
      }
      //   }
//   else{
//      int temp;
//      double q1, q2;
//      q1 = GetAlleleProbs( genotype[0], ancestry(0), AlleleFreqs ) * GetAlleleProbs( genotype[1], ancestry(1), AlleleFreqs );
//      q2 = GetAlleleProbs( genotype[1], ancestry(0), AlleleFreqs ) * GetAlleleProbs( genotype[1], ancestry(0) , AlleleFreqs);
//      hap(0) = genotype[0];
//      hap(1) = genotype[1];
//      if( myrand() > q1 / ( q1 + q2 ) ){
//         temp = hap(0);
//         hap(0) = hap(1);
//         hap(1) = temp;
//      }
//   }
   return( hap );
}

void CompositeLocus::SampleHaplotypePair(int hap[2],Vector_i Haplotypes, Vector_i ancestry)
{
  int i;
  double Probs[ 2*Haplotypes.GetNumberOfElements()];
  
   for( int k = 0; k < Haplotypes.GetNumberOfElements(); k++ ){
    Probs[ 2*k ]     = HaplotypeProbs( Haplotypes(k) )( ancestry(0), ancestry(1) );
    Probs[ 2*k + 1 ] = HaplotypeProbs( Haplotypes(k) )( ancestry(1), ancestry(0) );
  }
   i = SampleFromDiscrete3( Probs,  2*Haplotypes.GetNumberOfElements() );
  DipLoop.Reset();
  DipLoop.Increment( Haplotypes( i/2 ) );
  if( i % 2 ){
    hap[0] = (DipLoop.GetCountParent( 0 )*pmult).Sum();
    hap[1] = (DipLoop.GetCountParent( 1 )*pmult).Sum();
  }
  else{
    hap[1] = (DipLoop.GetCountParent( 0 )*pmult).Sum();
    hap[0] = (DipLoop.GetCountParent( 1 )*pmult).Sum();
  }
}

/**
 * Called every time the haplotype frequencies change. Constructs a 
 * three-dimensional matrix of Haplotype->Paternal Ancestry->Maternal
 * Ancestry, giving probabilities.
 */
void CompositeLocus::ConstructHaplotypeProbs(Matrix_d &AlleleFreqs)
{
   Vector_i DipBase, HapBase, x;
   Matrix_d GenoTypeProbs( Populations, Populations );

   DipBase.SetNumberOfElements( NumberOfLoci * 2 );
   HapBase.SetNumberOfElements( NumberOfLoci );
   pmult.SetNumberOfElements( NumberOfLoci );
   pmult.SetElements( 1 );
   x.SetNumberOfElements( NumberOfLoci * 2 );
   for( int k = 0; k < NumberOfLoci; k++ ){
      HapBase( k ) = NumberOfAlleles( k );
      DipBase( 2 * k ) = NumberOfAlleles( k );
      DipBase( 2 * k + 1 ) = NumberOfAlleles( k );
      for( int kk = 0; kk < k; kk++ ){
         pmult(kk) *= NumberOfAlleles( k );
      }
   }
   DipLoop.SetBase( DipBase );
   HapLoop.SetBase( HapBase );
   HaplotypeProbs.SetNumberOfElementsWithDimensions( DipLoop.GetDecimalBase(), Populations, Populations );

   float prob1, prob2;
   Vector_i count, mother, father;
   DipLoop.Reset();
   do{
      count = DipLoop.GetCount();
      father = DipLoop.GetCountParent( 0 );
      mother = DipLoop.GetCountParent( 1 );
      for( int pop1 = 0; pop1 < Populations; pop1++ ){
         for( int pop2 = 0; pop2 < Populations; pop2++ ){
            if( (father*pmult).Sum() < NumberOfStates - 1 )
               prob1 = AlleleFreqs( (father*pmult).Sum(), pop1 );
            else
               prob1 = 1 - (AlleleFreqs.GetColumn( pop1 )).Sum();
            if( (mother*pmult).Sum() < NumberOfStates - 1 )
               prob2 = AlleleFreqs( (mother*pmult).Sum(), pop2 );
            else
               prob2 = 1 - (AlleleFreqs.GetColumn( pop2 )).Sum();
            HaplotypeProbs( DipLoop.GetDecimalCount() )( pop1, pop2 ) = prob1 * prob2;
         }}
      DipLoop.Increment(1);
   }while( DipLoop.GetDecimalCount() != 0 );
}

/**
 * Given a list of possible haplotype pairs, returns sums of probabilities of these haplotypes
 * given each possible ordered pair of locus ancestry states 
 * 
 * Haplotypes - a list of possible haplotypes pairs compatible with the observed genotypes
 *
 * fixed - indicates whether the allelefrequencies are fixed
 * RandomAlleleFreqs - indicates whether the allele frequencies are random - why two indicators?
 *
 * returns:
 * a matrix with rows and cols indexing paternal and maternal ancestry. For 
 *   example, for African and European populations:
 *
 *       | AFR | EUR |
 *   ----|-----|-----|
 *   AFR | 0.5 | 0.2 |
 *   ----|-----|-----|
 *   EUR | 0.2 | 0.1 |
 *   ----|-----|-----|
 *
 *   n.b. the sum of all probabilities might not equal 1.0 - but
 *   probabilities are in correct proportions.
 */
//this matrix has to be rearranged as a vector for the HMM  
void CompositeLocus::GetGenotypeProbs(Matrix_d *Probs, Vector_i Haplotypes, bool fixed, int RandomAlleleFreqs)
{
  //Matrix_d GenoTypeProbs( Populations, Populations );
  Probs->SetNumberOfElements(Populations, Populations);// ! Possible Memory Leak !
   for( int k = 0; k < Haplotypes.GetNumberOfElements(); k++ ){
      if( fixed && RandomAlleleFreqs == 1 )
         *Probs += HaplotypeProbsMAP(Haplotypes(k) );
      else
         *Probs += HaplotypeProbs( Haplotypes(k) );
   }

   //return( GenoTypeProbs );
}

void CompositeLocus::GetGenotypeProbs(double **Probs, Vector_i Haplotypes, bool fixed, int RandomAlleleFreqs)
//Probs should be a k x k array to hold probabilities of genotype at this locus
{
  for(int i=0;i<Populations;++i)for(int j=0;j<Populations;++j)Probs[i][j]=0.0;//might be unnecessary if loops are rearranged
   for( int k = 0; k < Haplotypes.GetNumberOfElements(); k++ ){
      if( fixed && RandomAlleleFreqs == 1 )
	for(int i=0; i < Populations;++i)for(int j=0; j< Populations; ++j)
	  Probs[i][j] += HaplotypeProbsMAP(Haplotypes(k) )(i,j);
      else
	for(int i=0; i < Populations;++i)for(int j=0; j< Populations; ++j)
	  Probs[i][j] += HaplotypeProbs( Haplotypes(k) )(i,j);
   }

}

// doesn't really calculate posterior mode
// just sets to current value of hap freqs.  ok for Chib algorithm if strong prior
//this either misnamed or misdefined
void CompositeLocus::setHaplotypeProbsMAP()
{
   HaplotypeProbs = HaplotypeProbsMAP;
}

/**
 * Called only by method UpdateScoresForSNPsWithinHaplotype in ScoreTests
 * Given an unordered genotype, returns a Vector_i containing number of copies of allele 
 * 2 at each simple locus in the composite locus.
 * Used to test individual loci in haplotype for association.
 * 
 * genotype - a two-element STL vector in which each element is a one-dimensional array of
 * alleles coded as unsigned integers numbered starting at 0 ("decoded" format).  
 * a vector, of length equal to the number of simple loci in this composite
 * locus, containing the number of copies of allele 2 at each locus.
 *
 * n.b. this method is only useful in composite loci composed of diallelic simple loci
 * should be generalized to deal with multi-allelic loci
 */
Vector_i CompositeLocus::GetAlleleCountsInHaplotype(std::vector<unsigned short >&genotype)
{
  /**
   * AlleleCounts contains counts of the number of 2 alleles at each
   * locus in haplotype.  Used to test individual loci in haplotype
   * for association.  Only use for haplotypes made up of SNPs.
   */

   Vector_i AlleleCounts( NumberOfLoci );
   Vector_i decoded = decodeGenotype(genotype);// subtract 1 from allele numbers

   for( int k = 0; k < NumberOfLoci; k++ ){
      if(decoded(k*2)!=-1 && decoded(k*2+1)!=-1){
	if(decoded(k*2) == 1){
	  AlleleCounts(k)++;
	}
	if(decoded(k*2+1) == 1){
	  AlleleCounts(k)++;
	}
      } else {
	AlleleCounts(k) = 99;
      }
   }
   return( AlleleCounts );
}

// can get rid of this once we eliminate special methods for haploid data
// method takes a single encoded genotype as argument and subtracts 1 from the allele numbers 
// so that alleles are numbered from 0
// presumably missing genotypes will be recoded as pairs of minus ones 
Vector_i CompositeLocus::decodeGenotype(std::vector<unsigned short >&encoded)
{
  Vector_i decoded(encoded.size());

  for(unsigned int i=0;i<encoded.size();i++){
    decoded(i) = ((int)(encoded[i])) - 1;
  }
  return decoded;
}

Vector_i CompositeLocus::GetNumberOfAlleles()
{
  return NumberOfAlleles;
}

/**
 * Given an unordered genotype, calculates every possible haplotype pair that
 * is compatible with observed genotype.
 * Should be called once only for each (composite) locus and each individual. 
 * The lists of possible haplotypes are stored in the Individual objects to save recalculating each iteration.
 *
 * genotype - a vector of genotypes
 *   
 *
 * Haplotypes:
 * a list of possible haplotype pairs (in VectorLoop format).
 */
void CompositeLocus::SetPossibleHaplotypes(Vector_i *Haplotypes, std::vector<unsigned short >&genotype)
{
   int MissingValues;
   Vector_i LociBase, WhereMissingValues, count, xx, x;
   VectorLoop LociLoop;
   //Vector_i Haplotypes;

   x.SetNumberOfElements( 2*NumberOfLoci );
   LociBase.SetNumberOfElements( NumberOfLoci );
   WhereMissingValues.SetNumberOfElements( NumberOfLoci );

   x = decodeGenotype(genotype);
   
   MissingValues = 0;
   for( int k = 0; k < NumberOfLoci; k++ ){
      if( x(2*k) == -1 ){
         WhereMissingValues( MissingValues ) = k;
         MissingValues++;
      }
   }
   if( MissingValues == 0 ){
      Haplotypes->SetNumberOfElements( NumberOfStates );
      for( int jj = 0; jj < NumberOfStates; jj++ ){
         xx = Haplotype( x, jj );
         (*Haplotypes)(jj) = DipLoop.GetDecimal( xx );
      }
      Haplotypes->Distinct();
   }
   else{
      for( int k = 0; k < NumberOfLoci; k++ )
         if( x(2*k) != -1 )
            LociBase(k) = 2;
         else
            LociBase(k) = NumberOfAlleles(k) * NumberOfAlleles(k);
      
      LociLoop.SetBase( LociBase );
      Haplotypes->SetNumberOfElements( LociLoop.GetDecimalBase() );
      
      do{
         xx = x;
         count = LociLoop.GetCount();
         for( int k = 0; k < NumberOfLoci; k++ ){
            if( x(2*k) != -1 ){
               if( count(k) == 1 ){
                  xx(2*k) = x(2*k+1);
                  xx(2*k+1) = x(2*k);}}
            else{
               xx(2*k) = (int)(count(k)) / NumberOfAlleles(k);
               xx(2*k+1) = (int)count(k) % NumberOfAlleles(k);}}
	 (*Haplotypes)( LociLoop.GetDecimalCount() ) = DipLoop.GetDecimal( xx );
         LociLoop.Increment(1);
      }while( LociLoop.GetDecimalCount() != 0 );
      Haplotypes->Distinct();
   }
   //return Haplotypes;
}

int CompositeLocus::HapLoopGetDecimal(Vector_i x){
  return HapLoop.GetDecimal(x);
}

Vector_i CompositeLocus::Haplotype( Vector_i Genotype, int count )
{
   Vector_i hap;
   hap = Genotype;
   for( int i = 0; i < NumberOfLoci; i++ )
      if( (int)(count * pow(2.0, (double)(-i))) % 2 == 1 ){
         hap(2*i) = Genotype(2*i+1);
         hap(2*i+1) = Genotype(2*i);
      }
   
   return( hap );
}
// presumably this calculates score test for mis-spec allele freqs at multi-allelic loci
void CompositeLocus::UpdateScoreForMisSpecOfAlleleFreqs2(Matrix_d &AlleleFreqs, Matrix_i &AlleleCounts)
{
   double rn, r, pj, pi, q;
   Matrix_d NewScore( NumberOfStates - 1, 1 ), NewInfo( NumberOfStates - 1, NumberOfStates - 1 );
   for( int k = 0; k < Populations; k++ ){
      rn = (double)( AlleleCounts( NumberOfStates - 1, k ) );
      q = 1 - AlleleFreqs.GetColumn(k).Sum();
      for( int j = 0; j < NumberOfStates - 1; j++ ){
         r = AlleleCounts( j, k );
         pj = AlleleFreqs( j, k );
         NewScore( j, 0 ) = ( r / pj - rn / q ) * pj * ( 1 - pj );
         NewInfo( j, j ) = pj * ( 1 - pj )
            * ( r - ( rn / q ) * ( 2*pj - 1.0 - pj / q + pj * pj / q ) );
         for( int i = j+1; i < NumberOfStates - 1; i++ ){
            pi = AlleleFreqs( i, k );
            NewInfo( j, i ) = rn * pj * ( 1 - pj ) * pi * ( 1 - pi ) / ( q * q );
            NewInfo( i, j ) = NewInfo( j, i );
         }
      }
      SumNewScore(k) += NewScore;
      SumNewInfo(k) += NewInfo;
      SumNewScoreSq(k) += NewScore * NewScore.Transpose();
   }
}

/**
 * N.B. This only works for a single SNiP.
 * Updates what's required for the score tests. Only used with fixed
 * allele frequencies. This method is only used for monitoring.
 */
void CompositeLocus::UpdateScoreForMisSpecOfAlleleFreqs( Matrix_d phi, std::vector<unsigned short >&x, Matrix_d AlleleFreqs)
{
   Vector_d Pi( 3 ), Score( Populations );
   Pi.SetElements(0);
   for( int k = 0; k < Populations; k++ ){
      for( int kk = 0; kk < Populations; kk++ ){
         Pi(0) += AlleleFreqs( 0, k ) * AlleleFreqs( 0, kk ) * phi( k ,kk );
         Pi(1) += ( ( 1 - AlleleFreqs( 0, k ) ) * AlleleFreqs( 0, kk ) + ( 1 - AlleleFreqs( 0, kk ) ) * AlleleFreqs( 0, k ) ) * phi( k ,kk );
         Pi(2) += ( 1 - AlleleFreqs( 0, k ) ) * ( 1 - AlleleFreqs( 0, kk ) ) * phi( k ,kk );}}

   if( x[0] == 1 && x[1] == 1 ){
      for( int k = 0; k < Populations; k++ ){
         Score(k) = 2 * AlleleFreqs( 0, k ) * phi( k, k );
         for( int kk = 0; kk < Populations; kk++ )
            if( k != kk )
               Score(k) += AlleleFreqs( 0, kk ) * (phi( k, kk ) + phi( kk, k ));
         Score(k) /= Pi(0);
         ScoreGene( k, 0 ) += Score(k);
         InfoGene( k, k ) += Score(k) * Score(k) - 2 * phi( k, k ) / Pi(0);}
      for( int k = 0; k < Populations; k++ )
         for( int kk = 0; kk < Populations; kk++ )
             if( k != kk )
                InfoGene( k, kk ) += Score(k) * Score(kk) - (phi( k, kk ) + phi( kk, k )) / Pi(0);}
   
   else if( x[0] == 1 && x[1] != 1 ){
      for( int k = 0; k < Populations; k++ ){
         Score(k) = 2 * ( 1 - 2 * AlleleFreqs( 0, k ) ) * phi( k, k );
         for( int kk = 0; kk < Populations; kk++ )
            if( k != kk )
                Score(k) += ( 1 - 2 * AlleleFreqs( 0, kk ) ) * (phi( k, kk ) + phi( kk, k ));
         Score(k) /= Pi(1);
         ScoreGene( k, 0 ) += Score(k);
         InfoGene( k, k ) += Score(k) * Score(k) + 4 * phi( k, k ) / Pi(1);}
       for( int k = 0; k < Populations; k++ )
          for( int kk = 0; kk < Populations; kk++ )
             if( k != kk )
                InfoGene( k, kk ) += Score(k) * Score(kk) + 2*(phi( k, kk ) + phi( kk, k )) / Pi(1);}
   
   else if( x[0] != 0 && x[0] != 1 && x[1] != 1 ){
      for( int k = 0; k < Populations; k++ ){
          Score(k) = -2 * ( 1 - AlleleFreqs( 0, k ) ) * phi( k, k );
          for( int kk = 0; kk < Populations; kk++ )
             if( k != kk )
                Score(k) -= ( 1 - AlleleFreqs( 0, kk ) ) * (phi( k, kk ) + phi( kk, k ));
          Score(k) /= Pi(2);
          ScoreGene( k, 0 ) += Score(k);
          InfoGene( k, k ) += Score(k) * Score(k) - 2 * phi( k, k ) / Pi(2);}
      for( int k = 0; k < Populations; k++ )
         for( int kk = 0; kk < Populations; kk++ )
            if( k != kk )
               InfoGene( k, kk ) += Score(k) * Score(kk) - (phi( k, kk ) + phi( kk, k )) / Pi(2);}
}

/**
 * At the end of an iteration, zero's the score before calculating the
 * score afresh.  Only used with fixed allele frequencies.  This
 * method is only used for monitoring.  MisSpec eq 'mis-specification'
 */
void CompositeLocus::ResetScoreForMisSpecOfAlleleFreqs()
{
   ScoreGene.SetElements(0);
   InfoGene.SetElements(0);
}

/**
 * Keeps a running average of scores.
 * Only used with fixed allele frequencies. 
 * This method is only used for monitoring.
 */
void CompositeLocus::SumScoreForMisSpecOfAlleleFreqs()
{
   SumScoreGene += ScoreGene;
   SumInfoGene += InfoGene;
   SumScoreGeneSq += ScoreGene * ScoreGene.Transpose();
}

/**
 * 
 * Only used with fixed allele frequencies. 
 * This method is only used for monitoring.
 * TBC
 */
Matrix_d CompositeLocus::GetScore()
{
   return( SumScoreGene );
}

Matrix_d CompositeLocus::GetScoreSq()
{
   return( SumScoreGeneSq );
}

Matrix_d CompositeLocus::GetInfo()
{
   return( SumInfoGene );
}

Matrix_d CompositeLocus::GetNewScore( int k )
{
   return( SumNewScore(k) );
}

Matrix_d CompositeLocus::GetNewScoreSq( int k )
{
   return( SumNewScoreSq(k) );
}

Matrix_d CompositeLocus::GetNewInfo( int k )
{
   return( SumNewInfo(k) );
}

/**
 * Sets behaviour to not merge any haplotypes.
 * (Usually, rare haplotypes can be merged).
 * This method is only used for monitoring.
 */
void CompositeLocus::SetNoMergeHaplotypes()
{
   MergeHaplotypes.SetNumberOfElements( NumberOfStates );
   for( int i = 0; i < NumberOfStates; i++ )
      MergeHaplotypes(i) = i;
   NumberOfMergedHaplotypes = NumberOfStates;
}

/**
 * Decides which haplotypes to merge for score test, based on
 * frequencies <=0.01.  alpha is proportion of admix of each
 * population (used for weighting).  This method is only used for
 * monitoring.
 */
void CompositeLocus::SetDefaultMergeHaplotypes( Vector_d alpha, Matrix_d AlleleFreqs )
{
   int count = 0, count2 = 0;
   double p;
   Vector_i temp( NumberOfStates ), Merged( NumberOfStates );
   MergeHaplotypes.SetNumberOfElements( NumberOfStates );

   for( int i = 0; i < NumberOfStates - 1; i++ ){
      p = 0;
      for( int j = 0; j < Populations; j++ )
         p += alpha(j) * AlleleFreqs( i, j );
      p /= alpha.Sum();
      if( p > 0.01 ){
         temp( i ) = count;
         count++;
      }
      else
         temp(i) = NumberOfStates;
   }
   p = 0;
   for( int j = 0; j < Populations; j++ )
      p += alpha(j) * ( 1 - (AlleleFreqs.GetColumn( j )).Sum() );
   p /= alpha.Sum();
   if( p > 0.01 ){
      temp( NumberOfStates - 1 ) = count;
      count++;
   }
   else
      temp( NumberOfStates - 1 ) = NumberOfStates;

   for( int i = 0; i < NumberOfStates; i++ ){
      if( temp(i) == NumberOfStates )
         MergeHaplotypes(i) = count;
      else{
         MergeHaplotypes(i) = temp(i);
         Merged( count2 ) = i;
         count2++;
      }
   }

   if( count == NumberOfStates )
      NumberOfMergedHaplotypes = NumberOfStates;
   else
      NumberOfMergedHaplotypes = count + 1;

   VectorLoop hap( NumberOfAlleles );
   HapLabels.SetNumberOfElements( NumberOfMergedHaplotypes, NumberOfLoci );
   Vector_i temphaplabel;
   for( int i = 0; i < NumberOfMergedHaplotypes - 1; i++ ){
      hap.Reset();
      hap.Increment( Merged(i) );
      temphaplabel = hap.GetCount();
      HapLabels.SetRow(i, temphaplabel);
   }
   for( int j = 0; j < NumberOfLoci; j++ )
      HapLabels( NumberOfMergedHaplotypes - 1, j ) = 99;
}

Vector_i CompositeLocus::GetHapLabels( int i )
{
   return( HapLabels.GetRow(i) );
}

/**
 * Returns number of haplotypes that have been merged.
 * This method is only used for monitoring.
 */
int CompositeLocus::GetNumberOfMergedHaplotypes()
{
   return( NumberOfMergedHaplotypes );
}

/**
 * Given haplotype, returns merged haplotype.
 * This method is only used for monitoring.
 */
int CompositeLocus::GetMergedHaplotype( int i )
{
   return( MergeHaplotypes(i) );
}

// apparently calculates contribution of allele freqs to marginal likelihood of model
// by subtracting log prior density from log posterior
// in current version, this method is not called anywhere
double GetMarginalLikelihood( Vector_d PriorAlleleFreqs, Vector_d AlleleCounts )
{
   double f = gsl_sf_lngamma( PriorAlleleFreqs.Sum() ) -
      gsl_sf_lngamma( PriorAlleleFreqs.Sum() + AlleleCounts.Sum() );
   for( int i = 0; i < PriorAlleleFreqs.GetNumberOfElements(); i++ )
      f += gsl_sf_lngamma( PriorAlleleFreqs(i) + AlleleCounts(i) )
         - gsl_sf_lngamma( PriorAlleleFreqs(i) );
   return(f);
}






