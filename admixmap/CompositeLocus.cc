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
 * This is a ruler that measures 70 characters -
 * use it for manual line breaking.

         1         2         3         4         5         6         7
 2 4 6 8 0 2 4 6 8 0 2 4 6 8 0 2 4 6 8 0 2 4 6 8 0 2 4 6 8 0 2 4 6 8 0

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
  RandomAlleleFreqs = 0;
  Historical = false;
  
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
  AlleleFreqs = null_Matrix_d;
  HistoricalAlleleFreqs = null_Matrix_d;
  LikelihoodAlleleFreqs = null_Matrix_i;
  HistoricLikelihoodAlleleFreqs = null_Matrix_d;
  PriorAlleleFreqs = null_Matrix_d;
  SumAlleleFreqs = null_Matrix_d;
  SumEta = null_Vector_d;
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
  Fst = null_Vector_d;
  SumFst = null_Vector_d;
  SumNewScore = null_MatrixArray_d;
  SumNewInfo = null_MatrixArray_d;
  SumNewScoreSq = null_MatrixArray_d;
  HapLabels = null_Matrix_i;
}

CompositeLocus::~CompositeLocus()
{
   delete [] Label;
}

void
CompositeLocus::accept(LocusVisitor& v)
{
  v.visitCompositeLocus(*this);
}

/**
 * Changes number of loci in this composite locus, then sets each
 * locus to be diallelic (as for a SNP).
 *
 * NewNumberOfLoci - the number of loci to be represented by this
 * object.
 */
void
CompositeLocus::SetNumberOfLoci( int NewNumberOfLoci )
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
int
CompositeLocus::GetNumberOfLoci()
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
int
CompositeLocus::GetNumberOfStates()
{
   return( NumberOfStates );
}

/**
 * Indicates whether allele frequencies are fixed or random.
 *
 * Returns:
 * an integer representing a boolean. true (one) if allele frequencies
 * are random. false (zero) otherwise.
 */
int
CompositeLocus::IsRandom()
{
   return( RandomAlleleFreqs );
}

/**
 * Sets the number of alleles at a given locus.
 * Exits with an error if the locus doesn't exist.
 *
 * alleles - the number of alleles that exist at a given locus
 * locus - the given locus
 */
void
CompositeLocus::SetNumberOfAllelesOfLocus( int locus, int alleles )
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
int
CompositeLocus::GetNumberOfAllelesOfLocus( int locus )
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

void
CompositeLocus::AddLocus( int alleles )
{ 
   NumberOfAlleles.AddElement( NumberOfLoci );
   NumberOfAlleles( NumberOfLoci ) = alleles;
   NumberOfLoci++;
   NumberOfStates *= alleles;
}

void
CompositeLocus::SetNumberOfLabels()
{
   Label = new string[ NumberOfLoci ];
}

/**
 * Sets the frequencies of each allele at each locus in the composite
 * locus.
 *
 * NewAlleleFreqs - a two-dimensional matrix containing allele 
 *   frequencies. The first dimension is the allele number, being in
 *   the range of zero to two less than the number of states [see 
 *   GetNumberOfStates()]. The frequency of the final state is
 *   implied by the sum of all other frequencies. The second dimension
 *   is the population. Thus, for a composite locus with four states
 *   and European and African populations, the matrix might be:
 *
 *             Population
 *
 *            | EUR | AFR |
 *         ---|-----|-----|
 *          0 | 0.5 | 0.2 |
 *   State  1 | 0.1 | 0.2 |
 *          2 | 0.2 | 0.5 |
 *
 */

// should remove the Score and Info objects below from this class
void
CompositeLocus::SetAlleleFreqs( Matrix_d NewAlleleFreqs )
{
   if( NewAlleleFreqs.GetNumberOfRows() != NumberOfStates - 1 ){
      cout << "Error in number of alleles in SetAlleleFreqs.\n";
      cout << "Number of states = " << NumberOfStates << endl;
      cout << "AlleleFreqs has " << NewAlleleFreqs.GetNumberOfRows() << " rows.\n";
      exit(0);
   }

   AlleleFreqs = NewAlleleFreqs;
   Populations = AlleleFreqs.GetNumberOfCols();
   LikelihoodAlleleFreqs.SetNumberOfElements( NumberOfStates, Populations );

   ScoreGene.SetNumberOfElements( Populations, 1 );
   SumScoreGene.SetNumberOfElements( Populations, 1 );
   InfoGene.SetNumberOfElements( Populations, Populations );
   SumInfoGene.SetNumberOfElements( Populations, Populations );
   SumScoreGeneSq.SetNumberOfElements( Populations, Populations );

   SumNewScore.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, 1 );
   SumNewInfo.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, NumberOfStates - 1 );
   SumNewScoreSq.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, NumberOfStates - 1 );

   if( NumberOfLoci > 1 ){
      ConstructHaplotypeProbs();
      SetNoMergeHaplotypes();
   }
}

/**
 *
 * This method should be moved out of CompositeLocus into AlleleFreqs. 
 * Sets the frequencies of each allele at each locus in the
 * composite locus.
 *
 *   argument NewPriorAlleleFreqs - a two-dimensional matrix containing 
 *   parameters for the Dirichlet prior distribution of the allele frequencies. The first dimension is 
 *   the allele number, 
 *   being in the range of zero to two less than the number of states
 *   [see GetNumberOfStates()]. The sum of the prior parameters over all alleles in a population 
 *   (sumalpha) can be interpreted as 
 +   the "prior sample size". The second dimension is the population. Thus, for a 
 *   composite locus with four (? 5) states, and European and African 
 *   populations, the matrix might be:
 *
 *             Population
 *
 *            | EUR | AFR |
 *         ---|-----|-----|
 *          0 | 9.0 | 3.0 |
 *   State  1 | 3.0 | 4.0 |
 *          2 | 1.0 | 8.0 |
 *          3 | 2.0 | 1.0 |
 *
 *  AlleleFreqs is a matrix of expectations over this Dirichlet prior distribution, 
 *  calculated by dividing each prior parameter by the sum of the parameters. 
 * 
 */
void
CompositeLocus::SetPriorAlleleFreqs( Matrix_d NewPriorAlleleFreqs, bool fixed )
{
   double sumalpha;

   if( NewPriorAlleleFreqs.GetNumberOfRows() != NumberOfStates ){
      cout << "Error in number of alleles in SetPriorAlleleFreqs.\n";
      cout << "Number of states " << NumberOfStates << endl;
      cout << "PriorAlleleFreqs has " << NewPriorAlleleFreqs.GetNumberOfRows() << " rows.\n";
   }

   Populations = NewPriorAlleleFreqs.GetNumberOfCols();
   AlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );
   for( int j = 0; j < Populations; j++ ){
      sumalpha = ( NewPriorAlleleFreqs.GetColumn(j) ).Sum();
      for( int k = 0; k < NumberOfStates - 1; k++ )
         AlleleFreqs( k, j ) = ( NewPriorAlleleFreqs( k, j ) ) / sumalpha;
   }

   if( NumberOfLoci > 1 ){
      ConstructHaplotypeProbs();
      HaplotypeProbsMAP = HaplotypeProbs;
      SetNoMergeHaplotypes();
   }

   if( fixed ){
      ScoreGene.SetNumberOfElements( Populations, 1 );
      SumScoreGene.SetNumberOfElements( Populations, 1 );
      InfoGene.SetNumberOfElements( Populations, Populations );
      SumInfoGene.SetNumberOfElements( Populations, Populations );
      SumScoreGeneSq.SetNumberOfElements( Populations, Populations );
      SumNewScore.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, 1 );
      SumNewInfo.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, NumberOfStates - 1 );
      SumNewScoreSq.SetNumberOfElementsWithDimensions( Populations, NumberOfStates - 1, NumberOfStates - 1 );
      LikelihoodAlleleFreqs.SetNumberOfElements( NumberOfStates, Populations );
   }
   else{
      PriorAlleleFreqs = NewPriorAlleleFreqs;
      LikelihoodAlleleFreqs.SetNumberOfElements( NumberOfStates, Populations );
      SumAlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );
      RandomAlleleFreqs = 1;
   }
}

/**
 * This method also should be moved to AlleleFreqs class 
 * This method sets "historical allele frequencies", where the model has been specified to allow the 
 * allele freqs in the admixed population 
 * to vary from the historical allele frequencies in the unadmixed ancestral populations that have 
 * been sampled. 
 * Otherwise as for SetPriorAlleleFreqs
 * 
 */
void
CompositeLocus::SetHistoricalAlleleFreqs( Matrix_d NewHistoricLikelihoodAlleleFreqs )
{
   double sumalpha;

   if( NewHistoricLikelihoodAlleleFreqs.GetNumberOfRows() != NumberOfStates ){
      cout << "Error in number of alleles in SetHistoricalAlleleFreqs.\n";
      cout << "Number of states " << NumberOfStates << endl;
      cout << "HistoricalAlleleFreqs has "
           << NewHistoricLikelihoodAlleleFreqs.GetNumberOfRows() << " rows.\n";
   }
   Historical = true;
   HistoricLikelihoodAlleleFreqs = NewHistoricLikelihoodAlleleFreqs;
   PriorAlleleFreqs = HistoricLikelihoodAlleleFreqs + 0.501;
   Populations = HistoricLikelihoodAlleleFreqs.GetNumberOfCols();
   AlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );
   HistoricalAlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );
   LikelihoodAlleleFreqs.SetNumberOfElements( NumberOfStates, Populations );
   SumEta.SetNumberOfElements( Populations );

   for( int j = 0; j < Populations; j++ ){
      sumalpha = ( PriorAlleleFreqs.GetColumn(j) ).Sum();
      for( int k = 0; k < NumberOfStates - 1; k++ )
         AlleleFreqs( k, j ) = ( PriorAlleleFreqs( k, j ) ) / sumalpha;
   }

   if( NumberOfLoci > 1 ){
      ConstructHaplotypeProbs();
      SetNoMergeHaplotypes();
   }

   SumAlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );

   RandomAlleleFreqs = 1;
   Fst.SetNumberOfElements( Populations );
   SumFst.SetNumberOfElements( Populations );
   if( NumberOfStates > 2 ){
      MuProposal.resize( Populations );
      for( int k = 0; k < Populations; k++ ){
         MuProposal[k].SetParameters( 10, 0.01, 0.001, 0.1, 0.23 );
      }
   }
}

/**
 * Given the number of ancestral populations, sets default values for
 * allele frequencies and prior allele frequencies.
 *
 * populations - the number of ancestral populations
 */
void
CompositeLocus::SetDefaultAlleleFreqs( int populations )
{
   Populations = populations;

   if( Populations < 1 ){
      cout << "Error in SetDefaultAlleleFreqs( int populations ).\n";
      cout << "Number of populations = " << Populations << endl;
      exit(0);
   }
   if(NumberOfStates < 2){
     cout << "Error: The number of alleles at a locus is < 2. There must be at least two different alleles at each locus." << endl;
     exit(0);
   }

   PriorAlleleFreqs.SetNumberOfElements( NumberOfStates, Populations );
   PriorAlleleFreqs.SetElements( 0.5 );
   AlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );
   AlleleFreqs.SetElements( 1.0 / NumberOfStates );
   LikelihoodAlleleFreqs.SetNumberOfElements( NumberOfStates, Populations );

   if( NumberOfLoci > 1 ){
      ConstructHaplotypeProbs();
      SetNoMergeHaplotypes();
   }

   SumAlleleFreqs.SetNumberOfElements( NumberOfStates - 1, Populations );

   RandomAlleleFreqs = 1;
}

/**
 * Sets the name of this composite locus (usually from the 
 * allelefreqs.txt file
 *
 * newlabel - the name of this composite locus
 */
void
CompositeLocus::SetLabel( int index, string newlabel )
{
   Label[index] = newlabel;
}

/**
 * Gets the name of this composite locus
 *
 * returns:
 * the name of this composite locus
 */
string
CompositeLocus::GetLabel(int index)
{
   return( Label[index] );
}

/**
 * this method should be renamed UpdateAlleleCounts, as counts are sufficient stat but not likelihood
 * Given the unordered genotype and the ordered ancestry states at a
 * locus, this method randomly draws the phase of the genotype, then
 * updates the counts of alleles observed in each state of ancestry.
 * (MORE DETAILS PLEASE)
 */
void
CompositeLocus::UpdateLikelihoodAlleleFreqs(const vector<unsigned int>& genotype, Vector_i ancestry )
{
   Vector_i no, h(2);
   Vector_d Probs;
   Matrix_d ProbsM;

   if( NumberOfLoci == 1 ){
     if( genotype[0] ){ // no missing alleles
         ProbsM = GetLocusProbs(genotype,false);
         if( myrand() < ProbsM( ancestry(0), ancestry(1) ) / ( ProbsM( ancestry(0), ancestry(1) ) + ProbsM( ancestry(1), ancestry(0) ) ) ){
            LikelihoodAlleleFreqs( genotype[0] - 1, ancestry(0) )++;
            LikelihoodAlleleFreqs( genotype[1] - 1, ancestry(1) )++;
         }
         else{
            LikelihoodAlleleFreqs( genotype[0] - 1, ancestry(1) )++;
            LikelihoodAlleleFreqs( genotype[1] - 1, ancestry(0) )++;
         }
      }
   }

   else{
      h = SampleHaplotype( genotype, ancestry );
      LikelihoodAlleleFreqs( h(0), ancestry(1) )++;
      LikelihoodAlleleFreqs( h(1), ancestry(0) )++;
   }
}

void
CompositeLocus::UpdateLikelihoodAlleleFreqs_HaploidData(const vector<unsigned int>& genotype, int ancestry )
{
   int xx;
   if( NumberOfLoci == 1 )
      xx = genotype[0] - 1;
   else{
      Vector_i x = decodeGenotype(genotype);
      xx = HapLoop.GetDecimal( x );
   }
   LikelihoodAlleleFreqs( xx, ancestry )++;
}

/**
 * Given the genotype and ancestry of mother and father, returns a
 * haplotype pair.
 * 
 * genotype - a two-element vector of paternal and maternal genotypes
 *   in decimal notation (e.g. 1121 , 1221 ).
 * 
 * ancestry - a two-element vector of parental ancestry (e.g. 1,0 
 *   might represent european paternal and african maternal).
 *
 * returns:
 * two-element vector containing the haplotype pair
 */
Vector_i
CompositeLocus::SampleHaplotype(const vector<unsigned int>& genotype, Vector_i ancestry )
{
   Vector_i no, hap(2);
   if( NumberOfLoci > 1 ){
      int i;
      Vector_d Probs;
      no = GetHaplotypes( genotype );
      Probs.SetNumberOfElements( 2*no.GetNumberOfElements() );
      for( int k = 0; k < no.GetNumberOfElements(); k++ ){
         Probs( 2*k )     = HaplotypeProbs( no(k) )( ancestry(0), ancestry(1) );
         Probs( 2*k + 1 ) = HaplotypeProbs( no(k) )( ancestry(1), ancestry(0) );
      }
      i = SampleFromDiscrete2( &Probs );
      DipLoop.Reset();
      DipLoop.Increment( no( i/2 ) );
      if( i % 2 ){
         hap(0) = (DipLoop.GetCountParent( 0 )*pmult).Sum();
         hap(1) = (DipLoop.GetCountParent( 1 )*pmult).Sum();
      }
      else{
         hap(1) = (DipLoop.GetCountParent( 0 )*pmult).Sum();
         hap(0) = (DipLoop.GetCountParent( 1 )*pmult).Sum();
      }
   }
   else{
      int temp;
      double q1, q2;
      q1 = GetAlleleProbs( genotype[0], ancestry(0) ) * GetAlleleProbs( genotype[1], ancestry(1) );
      q2 = GetAlleleProbs( genotype[1], ancestry(0) ) * GetAlleleProbs( genotype[1], ancestry(0) );
      hap(0) = genotype[0];
      hap(1) = genotype[1];
      if( myrand() > q1 / ( q1 + q2 ) ){
         temp = hap(0);
         hap(0) = hap(1);
         hap(1) = temp;
      }
   }
   return( hap );
}

/**
 * LikelihoodAlleleFreqs is a matrix of counts of each allele in each
 * population this is the sufficient statistic for updating allele
 * frequencies, but shouldn't be called a likelihood.  Sets 
 * allele frequencies to 0 before they are updated.
 */
void
CompositeLocus::ResetLikelihoodAlleleFreqs()
{
   LikelihoodAlleleFreqs.SetElements(0);
}

Matrix_i
CompositeLocus::GetLikelihoodAlleleFreqs()
{
   return( LikelihoodAlleleFreqs );
}

Vector_i
CompositeLocus::GetLikelihoodAlleleFreqs( int population )
{
   return( LikelihoodAlleleFreqs.GetColumn( population ) );
}

void
CompositeLocus::SetLikelihoodAlleleFreqs( int, Matrix_i newLikelihoodAlleleFreqs )
{
   LikelihoodAlleleFreqs = newLikelihoodAlleleFreqs;
}

void
CompositeLocus::SetNumberOfStates( int newNumberOfStates )
{
   NumberOfStates= newNumberOfStates;
}

/**
 * This method should be moved to class AlleleFreqs
 * Whether the object should remember the results of sampling.
 *
 * flag - integer representing a boolean. Set true (one) to remember
 *   sampled data. Set false (zero) during burn in.
 */
void
CompositeLocus::SampleAlleleFreqs( int flag )
{
   Vector_d freqs;

   for( int j = 0; j < Populations; j++ ){
      freqs = gendirichlet( PriorAlleleFreqs.GetColumn(j)
                            + LikelihoodAlleleFreqs.GetColumn(j) );
      freqs.RemoveElement( NumberOfStates - 1 );
      AlleleFreqs.SetColumn( j, freqs );
      if( Historical ){
         freqs = gendirichlet( PriorAlleleFreqs.GetColumn(j)
                               + HistoricLikelihoodAlleleFreqs.GetColumn(j) );
         freqs.RemoveElement( NumberOfStates - 1 );
         HistoricalAlleleFreqs.SetColumn( j, freqs );
      }
   }

   if( NumberOfLoci > 1 )
      ConstructHaplotypeProbs();
   
   if( flag > 0 ){
      SumAlleleFreqs += AlleleFreqs;
   }
}

Matrix_d
CompositeLocus::GetSumAlleleFreqs()
{
   return( SumAlleleFreqs );
}

/**
 * Sets the sum of allele frequencies over all sampling iterations to
 * zero.
 */
void
CompositeLocus::ResetSumAlleleFreqs()
{
   SumAlleleFreqs.SetElements(0);
}

/**
 * Gets the frequencies of each haplotype in the composite locus.
 *
 * returns:
 * a matrix containing the frequencies of each allele
 * (what is the structure of this matrix - how many dimensions?)
 */
Matrix_d
CompositeLocus::GetAlleleFreqs()
{
   return( AlleleFreqs );
}

Vector_d
CompositeLocus::getAlleleFreqsMAP( int population )
{
   return( AlleleFreqsMAP.GetColumn( population ) );
}

/**
 * Returns Dirichlet parameters for allele frequencies for a particular population.
 * 
 * population - the number of the population (zero based)
 * 
 * returns:
 * a vector containing Dirichlet parameters for frequencies of each allele 
 * expected frequencies are calculated by dividing each parameter by the sum of parameters

 */
Vector_d
CompositeLocus::GetPriorAlleleFreqs( int population )
{
   return( PriorAlleleFreqs.GetColumn( population ) );
}

/**
 * Called every time the haplotype frequencies change. Constructs a 
 * three-dimensional matrix of Haplotype->Paternal Ancestry->Maternal
 * Ancestry, giving probabilities.
 */
void
CompositeLocus::ConstructHaplotypeProbs()
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

Matrix_d
CompositeLocus::GetLikelihood( const vector<unsigned int> genotype, bool diploid, bool fixed )
{
  Matrix_d Prob;
  if( diploid ){
     if( NumberOfLoci == 1 ){
        Prob = GetLocusProbs(genotype,fixed);
        if( genotype[0] != genotype[1] )
           for( int k = 0; k < Populations; k++ ){
              for( int kk = k; kk < Populations; kk++ ){
                 Prob(k,kk) = Prob(k,kk) + Prob(kk,k);
                 Prob(kk,k) = Prob(k,kk);
              }
           }
        
     } else {
        Prob = GetGenotypeProbs(genotype,fixed);
     }
  }
  else{
     Prob.SetNumberOfElements( Populations, 1 );
     if( NumberOfLoci == 1 ){
        for( int pop = 0; pop < Populations; pop++ ){
           Prob( pop, 0 ) = GetAlleleProbs( genotype[0] - 1, pop );
        }
     }
     else{
        Vector_i x = decodeGenotype(genotype);
        int xx = HapLoop.GetDecimal( x );
        for( int pop = 0; pop < Populations; pop++ ){
           Prob( pop, 0 ) = GetAlleleProbs( xx - 1, pop );
        }
     }
  }
  return( Prob );
}

/**
 * Given an unordered genotype, returns the probabilities of it having
 * an ancestry from a population.
 * 
 * genotype - a two-element vector of paternal and maternal genotypes.
 *   (how is this stored - in decimal or binary??)
 *
 * returns:
 * a two-dimensional matrix of paternal vs. maternal ancestry. For 
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
Matrix_d
CompositeLocus::GetGenotypeProbs(const vector<unsigned int>& genotype, bool fixed)
{
   Vector_i no;
   Matrix_d GenoTypeProbs( Populations, Populations );

   no = GetHaplotypes(genotype);
   for( int k = 0; k < no.GetNumberOfElements(); k++ ){
      if( fixed && RandomAlleleFreqs == 1 )
         GenoTypeProbs += HaplotypeProbsMAP( no(k) );
      else
         GenoTypeProbs += HaplotypeProbs( no(k) );
   }

   return( GenoTypeProbs );
}

void
CompositeLocus::setHaplotypeProbsMAP()
{
   HaplotypeProbs = HaplotypeProbsMAP;
}

/**
 * Given an unordered genotype, returns the number of times allele 
 * number-2 appears at each locus of the composite locus.
 * Used to test individual loci in haplotype for association.
 * 
 * genotype - a two-element vector of paternal and maternal genotypes
 *   in decimal notation (e.g. 1121 , 1221 ).
 *
 * returns:
 * a vector, of length equal to the number of loci in this composite
 * locus, containing the number of allele number-2 at each locus.
 *
 * n.b. this method is only useful in composite loci composed of SNPs
 */
Vector_i
CompositeLocus::GetAlleleCountsInHaplotype(const vector<unsigned int>& genotype)
{
  /**
   * AlleleCounts contains counts of the number of 2 alleles at each
   * locus in haplotype.  Used to test individual loci in haplotype
   * for association.  Only use for haplotypes made up of SNPs.
   */

   Vector_i AlleleCounts( NumberOfLoci );
   Vector_i decoded = decodeGenotype(genotype);

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

Vector_i
CompositeLocus::decodeGenotype(const vector<unsigned int>& encoded)
{
  Vector_i decoded(encoded.size());
  for(unsigned int i=0;i<encoded.size();i++){
    decoded(i) = ((int)(encoded[i])) - 1;
  }
  return decoded;
}

Vector_i
CompositeLocus::GetNumberOfAlleles()
{
  return NumberOfAlleles;
}

/**
 * Given an unordered genotype, returns every possible haplotype that
 * is possible.
 *
 * genotype - a two-element vector of paternal and maternal genotypes
 *   in decimal notation (e.g. 1121 , 1221 ).
 *
 * returns:
 * a list of possible haplotypes (in VectorLoop format).
 */
Vector_i
CompositeLocus::GetHaplotypes(const vector<unsigned int>& genotype)
{
   int MissingValues;
   Vector_i LociBase, no, WhereMissingValues, count, xx, x;
   VectorLoop LociLoop;

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
      no.SetNumberOfElements( NumberOfStates );
      for( int jj = 0; jj < NumberOfStates; jj++ ){
         xx = Haplotype( x, jj );
         no(jj) = DipLoop.GetDecimal( xx );
      }
      no.Distinct();
   }
   else{
      for( int k = 0; k < NumberOfLoci; k++ )
         if( x(2*k) != -1 )
            LociBase(k) = 2;
         else
            LociBase(k) = NumberOfAlleles(k) * NumberOfAlleles(k);
      
      LociLoop.SetBase( LociBase );
      no.SetNumberOfElements( LociLoop.GetDecimalBase() );
      
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
            no( LociLoop.GetDecimalCount() ) = DipLoop.GetDecimal( xx );
         LociLoop.Increment(1);
      }while( LociLoop.GetDecimalCount() != 0 );
      no.Distinct();
   }

   return( no );
}
//new version of GetHaplotypes to be called once only to set PossHaplotypes
void CompositeLocus::GetPossibleHaplotypes(vector< unsigned int >& genotype)
{
   int MissingValues;
   Vector_i LociBase,  WhereMissingValues, count, xx, x;
   VectorLoop LociLoop;

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
      PossHaplotypes.SetNumberOfElements( NumberOfStates );
      for( int jj = 0; jj < NumberOfStates; jj++ ){
         xx = Haplotype( x, jj );
         PossHaplotypes(jj) = DipLoop.GetDecimal( xx );
      }
      PossHaplotypes.Distinct();
   }
   else{
      for( int k = 0; k < NumberOfLoci; k++ )
         if( x(2*k) != -1 )
            LociBase(k) = 2;
         else
            LociBase(k) = NumberOfAlleles(k) * NumberOfAlleles(k);
      
      LociLoop.SetBase( LociBase );
      PossHaplotypes.SetNumberOfElements( LociLoop.GetDecimalBase() );
      
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
            PossHaplotypes( LociLoop.GetDecimalCount() ) = DipLoop.GetDecimal( xx );
         LociLoop.Increment(1);
      }while( LociLoop.GetDecimalCount() != 0 );
      PossHaplotypes.Distinct();
   }

}

/**
 * N.B. This only works with a simple locus.
 * Given an unordered genotype, returns a matrix representing the
 * probability of locus ancestry.
 */
Matrix_d
CompositeLocus::GetLocusProbs(const vector<unsigned int>& x, bool fixed )
{
   MatrixArray_d Prob( 2, Populations, 1 );
   
   for( int pop = 0; pop < Populations; pop++ )
   {
      for( int i = 0; i < 2; i++ )
      {
         if( fixed && RandomAlleleFreqs == 1 )
            Prob(i)( pop, 0 ) = GetAlleleProbsMAP( x[i]-1, pop );
         else
            Prob(i)( pop, 0 ) = GetAlleleProbs( x[i]-1, pop );
      }
   }

   return Prob(0) * Prob(1).Transpose();

}

double
CompositeLocus::GetAlleleProbs( int x, int ancestry )
{
   double P;
   if( x < NumberOfAlleles(0) - 1 )
      P = AlleleFreqs( x, ancestry );
   else
   {
      P = 1;
      for( int j = 0; j < NumberOfAlleles(0) - 1; j++ )
         P -= AlleleFreqs( j, ancestry );
   }
   return P;
}

double
CompositeLocus::GetAlleleProbsMAP( int x, int ancestry )
{
   double P;
   if( x < NumberOfAlleles(0) - 1 )
      P = AlleleFreqsMAP( x, ancestry );
   else
   {
      P = 1;
      for( int j = 0; j < NumberOfAlleles(0) - 1; j++ )
         P -= AlleleFreqsMAP( j, ancestry );
   }
   return P;
}

void
CompositeLocus::setAlleleFreqsMAP()
{
   AlleleFreqsMAP = AlleleFreqs;
}

Vector_i
CompositeLocus::Haplotype( Vector_i Genotype, int count )
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

void
CompositeLocus::UpdateScoreForMisSpecOfAlleleFreqs2()
{
   double rn, r, pj, pi, q;
   Matrix_d NewScore( NumberOfStates - 1, 1 ), NewInfo( NumberOfStates - 1, NumberOfStates - 1 );
   for( int k = 0; k < Populations; k++ ){
      rn = (double)( LikelihoodAlleleFreqs( NumberOfStates - 1, k ) );
      q = 1 - AlleleFreqs.GetColumn(k).Sum();
      for( int j = 0; j < NumberOfStates - 1; j++ ){
         r = LikelihoodAlleleFreqs( j, k );
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
void
CompositeLocus::UpdateScoreForMisSpecOfAlleleFreqs( Matrix_d phi, vector<unsigned int> x)
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
void
CompositeLocus::ResetScoreForMisSpecOfAlleleFreqs()
{
   ScoreGene.SetElements(0);
   InfoGene.SetElements(0);
}

/**
 * Keeps a running average of scores.
 * Only used with fixed allele frequencies. 
 * This method is only used for monitoring.
 */
void
CompositeLocus::SumScoreForMisSpecOfAlleleFreqs()
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
Matrix_d
CompositeLocus::GetScore()
{
   return( SumScoreGene );
}

Matrix_d
CompositeLocus::GetScoreSq()
{
   return( SumScoreGeneSq );
}

Matrix_d
CompositeLocus::GetInfo()
{
   return( SumInfoGene );
}

Matrix_d
CompositeLocus::GetNewScore( int k )
{
   return( SumNewScore(k) );
}

Matrix_d
CompositeLocus::GetNewScoreSq( int k )
{
   return( SumNewScoreSq(k) );
}

Matrix_d
CompositeLocus::GetNewInfo( int k )
{
   return( SumNewInfo(k) );
}

int
CompositeLocus::GetSize()
{
  return 1;
}

/**
 * Sets behaviour to not merge any haplotypes.
 * (Usually, rare haplotypes can be merged).
 * This method is only used for monitoring.
 */
void
CompositeLocus::SetNoMergeHaplotypes()
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
void
CompositeLocus::SetDefaultMergeHaplotypes( Vector_d alpha )
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

Vector_i
CompositeLocus::GetHapLabels( int i )
{
   return( HapLabels.GetRow(i) );
}

/**
 * Returns number of haplotypes that have been merged.
 * This method is only used for monitoring.
 */
int
CompositeLocus::GetNumberOfMergedHaplotypes()
{
   return( NumberOfMergedHaplotypes );
}

/**
 * Given haplotype, returns merged haplotype.
 * This method is only used for monitoring.
 */
int
CompositeLocus::GetMergedHaplotype( int i )
{
   return( MergeHaplotypes(i) );
}

/**
 * Used when model is specified to allow allele freqs in admixed
 * population to vary from the "historical" allele freqs in the
 * unadmixed population.  Dirichlet distribution for allele freqs at
 * locus with k alleles is specified with (k - 1) frequency parameters
 * (mu) and with a single dispersion parameter (eta) this method
 * samples mu and eta, and updates Dirichlet parameters for the allele
 * frequencies
 */
Vector_d
CompositeLocus::GetStatsForEta( int population )
{
   Vector_d stats( NumberOfStates );
   for( int i = 0; i < NumberOfStates - 1; i++ )
      stats( i ) = log( AlleleFreqs( i, population ) ) + log( HistoricalAlleleFreqs( i, population ) );
   stats( NumberOfStates - 1 ) = log( 1 - AlleleFreqs.GetColumn( population ).Sum() )
      + log( 1 - HistoricalAlleleFreqs.GetColumn( population ).Sum() );
   return stats;
}

void
CompositeLocus::UpdatePriorAlleleFreqs( int j, const Vector_d& mu )
{
   PriorAlleleFreqs.SetColumn( j, mu );
//    double sum;
//    Vector_d freqs;
//    for( int j = 0; j < Populations; j++ ){
//       freqs = PriorAlleleFreqs.GetColumn(j);
//       sum = freqs.Sum();
//       PriorAlleleFreqs.SetColumn( j, freqs * eta(j) / sum );
//    }
}

// Method samples the prior frequency parameters mu.
// Takes eta, the sum of the frequency parameters for each locus.
void
CompositeLocus::SamplePriorAlleleFreqs( Vector_d eta )
{
   if( NumberOfStates == 2 )
      SamplePriorAlleleFreqs1D( eta );
  else
     SamplePriorAlleleFreqsMultiDim( eta );
}

void
CompositeLocus::SamplePriorAlleleFreqsMultiDim( Vector_d eta )
{
   vector<int> accept(Populations,0);
   Vector_d mu1, mu2;
   for( int j = 0; j < Populations; j++ ){
      double Proposal1=0, Proposal2=0, f1=0, f2=0;
      mu1 = PriorAlleleFreqs.GetColumn(j) / eta(j);
      mu2 = gendirichlet( mu1 / MuProposal[j].GetSigma() );
      
      for( int i = 0; i < NumberOfStates; i++ ){
         f1 += 0.1 * log( mu1(i) ) + 0.1 * log( 1 - mu1(i) );
         f2 += 0.1 * log( mu2(i) ) + 0.1 * log( 1 - mu2(i) );
         Proposal1 += (eta(j) * mu2(i) - 1) * log( mu1(i) ) - gsl_sf_lngamma( eta(j) * mu2(i) );
         Proposal2 += (eta(j) * mu1(i) - 1) * log( mu2(i) ) - gsl_sf_lngamma( eta(j) * mu1(i) );
      }
      
      int numberofstates = mu1.GetNumberOfElements();
      for( int k = 0; k < numberofstates; k++ ){
         f1 -= Populations * gsl_sf_lngamma( mu1( k )* eta(j) );
         f2 -= Populations * gsl_sf_lngamma( mu2( k )* eta(j) );
            f1 += gsl_sf_lngamma( mu1( k )* eta(j) + LikelihoodAlleleFreqs(k,j) );
            f2 += gsl_sf_lngamma( mu2( k )* eta(j) + LikelihoodAlleleFreqs(k,j) );
            f1 += gsl_sf_lngamma( mu1( k )* eta(j) + HistoricLikelihoodAlleleFreqs(k,j) );
            f2 += gsl_sf_lngamma( mu2( k )* eta(j) + HistoricLikelihoodAlleleFreqs(k,j) );
      }
      if( log(myrand()) < f2 - f1 - Proposal2 + Proposal1 ){
         PriorAlleleFreqs.SetColumn( j, mu2 * eta(j) );
         accept[j] = 1;
         MuProposal[j].Event(true);
      }
      else
         MuProposal[j].Event(false);
   }
}

void
CompositeLocus::SamplePriorAlleleFreqs1D( Vector_d eta )
{
   double lefttruncation = 0.1;
   Vector_d MuParameters(2);
   MatrixArray_i counts0(1);
   MatrixArray_d counts1(1);

// Construct adaptive rejection sampler for mu.
   counts0(0) = LikelihoodAlleleFreqs;
   counts1(0) = HistoricLikelihoodAlleleFreqs;
   //warning message - move to wherever loci data are read
   for(int row=0;row<counts0(0).GetNumberOfRows();++row)for(int col=0;col<counts0(0).GetNumberOfCols();++col)
     if(counts0(0)(row,col)==0 && counts1(0)(row,col)==0){
       //       cout<<"Warning: zero copies of allele ("<<row<<","<<col<<") in both admixed and unadmixed samples"<<endl;
   //poss return name of comp locus
     }
   DARS SampleMu( 0, 0, 0, MuParameters, fMu, dfMu, ddfMu, counts0, counts1 );

   SampleMu.SetLeftTruncation( lefttruncation );
   for( int j = 0; j < Populations; j++ ){
//      sum = 0.0;
      MuParameters(0) = eta(j);
      MuParameters(1) = j;
//       MuParameters(1) = (float)LikelihoodAlleleFreqs( NumberOfStates - 1, j );
//       MuParameters(2) = (float)HistoricLikelihoodAlleleFreqs( NumberOfStates - 1, j );
//       MuParameters(3) = (float)LikelihoodAlleleFreqs( k, j );
//       MuParameters(4) = (float)HistoricLikelihoodAlleleFreqs( k, j );
      SampleMu.SetRightTruncation( eta(j) - lefttruncation );
      SampleMu.UpdateParameters( MuParameters );
      PriorAlleleFreqs( 0, j ) = SampleMu.Sample();
// Last prior frequency parameter is determined by; sum of mu's = eta.
      PriorAlleleFreqs( 1, j ) = eta(j) - PriorAlleleFreqs( 0, j );
   }
}

double
GetMarginalLikelihood( Vector_d PriorAlleleFreqs, Vector_d LikelihoodAlleleFreqs )
{
   double f = gsl_sf_lngamma( PriorAlleleFreqs.Sum() ) -
      gsl_sf_lngamma( PriorAlleleFreqs.Sum() + LikelihoodAlleleFreqs.Sum() );
   for( int i = 0; i < PriorAlleleFreqs.GetNumberOfElements(); i++ )
      f += gsl_sf_lngamma( PriorAlleleFreqs(i) + LikelihoodAlleleFreqs(i) )
         - gsl_sf_lngamma( PriorAlleleFreqs(i) );
   return(f);
}

double
fMu( Vector_d &parameters, MatrixArray_i& counts0, MatrixArray_d& counts1, double mu )
{
//    Vector_d rn(2), ri(2);
//    rn(0) = parameters(1);
//    rn(1) = parameters(2);
//    ri(0) = parameters(3);
//    ri(1) = parameters(4);
   int pop = (int)parameters(1);
   double eta = parameters(0);
   double prior = 0.1 * log( mu / eta ) + 0.1 * log( 1 - mu / eta );
   double f = prior - 2 * gsl_sf_lngamma( mu ) - 2 * gsl_sf_lngamma( eta - mu );
//   for( int i = 0; i < 2; i++ )
//      f += gsl_sf_lngamma( mu + ri(i) ) + gsl_sf_lngamma( eta - mu + rn(i) );
   f += gsl_sf_lngamma( mu+counts0(0)(0,pop) ) + gsl_sf_lngamma( eta-mu+counts0(0)(1,pop) );
   f += gsl_sf_lngamma( mu+counts1(0)(0,pop) ) + gsl_sf_lngamma( eta-mu+counts1(0)(1,pop) );

   return f;
}

double
dfMu( Vector_d &parameters, MatrixArray_i& counts0, MatrixArray_d& counts1, double mu )
{
//    Vector_d rn(2), ri(2);
//    rn(0) = parameters(1);
//    rn(1) = parameters(2);
//    ri(0) = parameters(3);
//    ri(1) = parameters(4);
   int pop = (int)parameters(1);
   double eta = parameters(0), x, y1, y2;
   double prior = 0.1 / mu - 0.1 / ( eta - mu );
   double f = prior;
   x = parameters(0) - mu;
   if(mu < 0)cout<<"\nError in dfMu in compositelocus.cc - arg mu to ddigam is negative\n"; 
   ddigam( &mu, &y1 );
   if(x < 0)cout<<"\nError in dfMu in compositelocus.cc - arg x to ddigam is negative\n"; 
   ddigam( &x, &y2 );
   f += 2 * ( y2 - y1 );

//    for( int i = 0; i < 2; i++ ){
//       x = mu + ri(i);
//       ddigam( &x, &y2 );
//       f += y2;
//       x = eta - mu + rn(i);
//       ddigam( &x, &y2 );
//       f -= y2;
//    }

   x = mu + counts0(0)(0,pop);
   ddigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts0(0)(1,pop);
   ddigam( &x, &y2 );
   f -= y2;

   x = mu + counts1(0)(0,pop);
   ddigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts1(0)(1,pop);
   ddigam( &x, &y2 );
   f -= y2;

   return f;
}

double
ddfMu( Vector_d &parameters, MatrixArray_i& counts0, MatrixArray_d& counts1, double mu )
{
//    Vector_d rn(2), ri(2);
//    rn(0) = parameters(1);
//    rn(1) = parameters(2);
//    ri(0) = parameters(3);
//    ri(1) = parameters(4);
   int pop = (int)parameters(1);
   double eta = parameters(0), x, y1, y2;
   double prior = -0.1 / (mu*mu) - 0.1 / (( eta - mu ) * ( eta - mu ) );
   double f = prior;
   x = parameters(0) - mu;
   trigam( &mu, &y1 );
   trigam( &x, &y2 );
   f -= 2 * ( y2 + y1 );

//    for( int i = 0; i < 2; i++ ){
//       x = mu + ri(i);
//       trigam( &x, &y2 );
//       f += y2;
//       x = eta - mu + rn(i);
//       trigam( &x, &y2 );
//       f += y2;
//    }

   x = mu + counts0(0)(0,pop);
   trigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts0(0)(1,pop);
   trigam( &x, &y2 );
   f += y2;

   x = mu + counts1(0)(0,pop);
   trigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts1(0)(1,pop);
   trigam( &x, &y2 );
   f += y2;

   return f;
}

void
CompositeLocus::UpdateFst()
{
   double q_admix,q_parental,f,H_admix, H_parental, H_combined, pbar;
   for( int k = 0; k < Populations; k++ ){
      H_admix = 0;
      H_parental = 0;
      H_combined = 0;

      for( int i = 0; i < NumberOfStates - 1; i++ ){
         H_admix += AlleleFreqs( i, k ) * AlleleFreqs( i, k );
         H_parental += HistoricalAlleleFreqs( i, k ) * HistoricalAlleleFreqs( i, k );
         pbar = 0.5 * ( AlleleFreqs( i, k ) + HistoricalAlleleFreqs( i, k ) );
         H_combined += pbar * pbar;
      }

      q_admix = 1 - AlleleFreqs.GetColumn(k).Sum();
      H_admix += q_admix * q_admix;
      q_parental = 1 - HistoricalAlleleFreqs.GetColumn(k).Sum();
      H_parental += q_parental * q_parental;
      pbar = 0.5 * ( q_admix + q_parental );
      H_combined += pbar * pbar;

      H_combined = 1 - H_combined;
      H_admix = 1 - H_admix;
      H_parental = 1 - H_parental;
      f = ( H_combined - 0.5 * ( H_admix + H_parental ) ) / H_combined;
      Fst( k ) = 2*f / ( 1 + f );
   }
   SumFst += Fst;
}

Vector_d
CompositeLocus::GetFst()
{
   return( Fst );
}

