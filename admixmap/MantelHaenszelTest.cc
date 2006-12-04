#include "MantelHaenszelTest.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "Comms.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_math.h"//for gsl_finite

MantelHaenszelTest::MantelHaenszelTest(){

}
MantelHaenszelTest::~MantelHaenszelTest(){

}

void MantelHaenszelTest::Initialise(unsigned NumStates, unsigned NumLoci){
  K = NumStates;
  Ksq = K*K;

  //assign K^2 2x2 tables
  //NOTE: asssuming all SNPs
  for(unsigned k = 0; k < Ksq; ++k){
    std::vector<unsigned> v(4);
    CountTable.push_back(v);
  }
  for(unsigned i = 0; i < NumLoci; ++i){
    Score.push_back(0);
    ScoreSq.push_back(0);
    Info.push_back(0);
  }
}

void MantelHaenszelTest::Update(const IndividualCollection* IC, const Genome& Loci){
  const std::vector<unsigned> nulltable(4, 0);
  const int NumWorkers = Comms::getNumWorkers();
  int ancA[2];//ancestry at A
  int ancB[2];//ancestry at B
  int locus = 0;

  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    for(unsigned j = 0; j < Loci.GetSizeOfChromosome(c)-1; ++j){
      //reset counts to zero
      fill(CountTable.begin(), CountTable.end(), nulltable);
      std::vector<unsigned> N(Ksq, 0);//number of gametes
      for(int i = Comms::getWorkerRank(); i < IC->getSize(); i += NumWorkers) {
	Individual* ind = IC->getIndividual(i);
	if( !ind->GenotypeIsMissing(locus) && !ind->GenotypeIsMissing(locus+1) ) {
	  //skip missing genotypes as hap pairs not sampled
	  ind->GetLocusAncestry(c, j, ancA);
	  ind->GetLocusAncestry(c, j+1, ancB);
	  const int* hA = ind->getSampledHapPair(locus);//realized hap pair at locus A
	  const int* hB = ind->getSampledHapPair(locus+1);//realized hap pair at locus B
	  int numGametes = 2;//number of gametes per 'individual'
	  if(hA[1] < 0) numGametes = 1;//haplotypes are coded as happair with second hap as -1
	  for(int g = 0; g < numGametes; ++g) {
	    ++CountTable[ancA[g]*K + ancB[g]][hA[g]*2 + hB[g]];
	    ++N[ancA[g]*K + ancB[g]];
	  }//end gamete loop
	}
      }//end individual loop

      for(unsigned k = 0; k < Ksq; ++k){
	if( N[k] > 0){//skip empty tables
	  const unsigned rowsum0 = CountTable[k][0]+CountTable[k][1];//sum of 1st row
	  if(rowsum0 ==0) continue;
	  const unsigned rowsum1 = CountTable[k][2]+CountTable[k][3];//       2nd
	  if(rowsum1 ==0) continue;
	  const unsigned colsum0 = CountTable[k][0]+CountTable[k][2];//sum of 1st col
	  if(colsum0 ==0) continue;
	  const unsigned colsum1 = CountTable[k][1]+CountTable[k][3];//       2nd
	  if(colsum1 ==0) continue;

	  const double score = (double)CountTable[k][0] - (double)( rowsum0*colsum0 ) / (double) N[k];
	  Score[locus] += score;
	  ScoreSq[locus] += score*score;
	  Info[locus] += (double)(rowsum0*rowsum1*colsum0*colsum1) / (double) (N[k]*N[k] *(N[k]-1));
	}
      }
      ++locus;
    }
    ++locus;//skip last locus on chromosome
  }//end locus loop
  
}

void MantelHaenszelTest::Output(const char* filename, const Genome& Loci, unsigned NumIters, const std::vector<std::string>& LocusLabels){
  outfile.open(filename);
  //write header
  outfile << "Loci\tScore\tCompInfo\tObsInfo\tPercentInfo\tzscore\tPValue\n";

  unsigned locus = 0;
  //loop over pairs of loci
  //NOTE: Exp = sum_over_tables ((obs first cell) - (Exp first cell))
  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    for(unsigned j = 0; j < Loci.GetSizeOfChromosome(c)-1; ++j){
      std::string label = LocusLabels[locus] + "/" + LocusLabels[locus+1];
      OutputScalarScoreTest(NumIters, &outfile, label, Score[locus], ScoreSq[locus], Info[locus], true);


//       const double ebar = Score[locus] / (double)NumIters;
//       const double vbar = Info[locus] / (double)NumIters;
//       double chisq = ( ebar)*(ebar) / vbar;
//       outfile << "\"" << LocusLabels[locus] << "/" << LocusLabels[locus+1] << "\"\t" 
// 	      << ebar << "\t" << vbar << "\t" << chisq << "\t";
//       if(gsl_finite(chisq)){
// 	double pvalue = gsl_cdf_chisq_Q (chisq, 1);
// 	outfile << pvalue << endl;
//       }
//       else
// 	outfile << "NA" << endl;
       ++locus;
    }
    ++locus;//skip last locus on chromosome
  }

  outfile.close();
}
