/*
  Class implements a score test for deviation from Hardy-Weinberg equilibrium
  in order to test for genotyping errors.
  This version evaluates for each single locus and sums over individuals.
  code for evaluation for each individual (rather than summing) is commented out.
*/
#include "HWTest.h"
#include "functions.h"
#include "gsl/gsl_cdf.h"

HWTest::HWTest(){
  score = 0;
  sumscore = 0;
  sumscore2 = 0;
  suminfo = 0;
  NumInd = 0;
  NumLoci = 0;
  samples  = 0;
}

HWTest::~HWTest(){
  //free_matrix(sumscore, NumInd);
  //free_matrix(sumscore2, NumInd);
  //free_matrix(suminfo, NumInd);
  delete[] score;
  delete[] sumscore;
  delete[] sumscore2;
  delete[] suminfo;
}

//void HWTest::Initialise(AdmixOptions *options, int nind, int nloci, LogWriter *Log){
void HWTest::Initialise(const AdmixOptions* const options, int nloci, LogWriter &Log){

  //NumInd = nind;
  NumLoci = nloci;
  Log.setDisplayMode(Quiet);

  if( options->getHWTestIndicator() ){//not really necessary
    if ( strlen( options->getHWTestFilename() ) ){
      outputfile.open( options->getHWTestFilename(), ios::out );
      if( !outputfile ){
	Log.setDisplayMode(Quiet);
	Log << "ERROR: Couldn't open hwtestfile\n";
	exit( 1 );}
      else {
	Log << "HW test file: " << options->getHWTestFilename() << "\n";
	
	//sumscore = alloc2D_d(NumInd, NumLoci);
	//sumscore2 = alloc2D_d(NumInd, NumLoci);
	//suminfo = alloc2D_d(NumInd, NumLoci);
	score = new double[NumLoci];
	sumscore = new double[NumLoci];
	sumscore2 = new double[NumLoci];
	suminfo = new double[NumLoci];
	for(int i = 0; i < NumLoci; ++i){
	  sumscore[i] = 0.0;
	  sumscore2[i] = 0.0;
	  suminfo[i] = 0.0;
	}
      }
    }
    else{
      Log << "No hwtestfile given\n";
      //exit(1);}
    }
  }
}    

void HWTest::Update(const IndividualCollection* const IC, const Chromosome* const* C, const Genome* const Loci){
  double H; // prob heterozygous
  bool h;
  int locus = 0, complocus = 0, Ancestry0, Ancestry1;
  vector<vector<unsigned short> > genotype;
  Individual *ind = 0;
  double **Prob0 = 0, **Prob1 = 0;

  //reset score
  for(int j = 0; j < NumLoci; ++j){
    score[j] = 0.0;
  }
  for(int i = 0; i < IC->getSize(); ++i){

    ind = IC->getIndividual(i);
    locus = 0;    //absolute locus number
    complocus = 0;//absolute comp locus number
    for(unsigned int chr = 0; chr < Loci->GetNumberOfChromosomes(); ++chr){ //loop over chromosomes
      for(unsigned int j = 0; j < C[chr]->GetSize(); ++j){                  //loop over comp loci on chromosome
	Ancestry0 = ind->GetLocusAncestry( (int)chr, 0, j);
	Ancestry1 = ind->GetLocusAncestry( (int)chr, 1, j);
	// if(Ancestry0 == Ancestry1) {

	  genotype = ind->getGenotype(complocus);

	  //allocate arrays to hold marginal alleleprobs; could be done in function but easier to control here.      
	  Prob0 = new double*[(*Loci)(complocus)->GetNumberOfLoci()];
	  Prob1 = new double*[(*Loci)(complocus)->GetNumberOfLoci()];
	  
	  for(int jj = 0; jj < (*Loci)(complocus)->GetNumberOfLoci(); ++jj){
	    Prob0[jj] = new double[ (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj)];
	    Prob1[jj] = new double[ (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj)];
	    
	  }
	  //retrieve marginal alleleprobs for composite locus complocus, given current ancestry states on each gamete of ind
	  (*Loci)(complocus)->getLocusAlleleProbs(Prob0, Ancestry0);// #loci x #alleles array  
	  (*Loci)(complocus)->getLocusAlleleProbs(Prob1, Ancestry1);   
	  
	  for(int jj = 0; jj < Loci->getNumberOfLoci(complocus); ++jj){       //loop over loci within comp locus
	    if( genotype[jj][0] != 0){ //non-missing genotype, assumes second gamete missing if first is
	      h = (genotype[jj][0] != genotype[jj][1]);
	      H = 1.0;
	      //compute prob of heterozygosity by subtracting from 1 the prob of homozygosity, ie sum of diagonal products	    
	      for(int a = 0; a < (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj); ++a){//loop over alleles
		H -= Prob0[jj][a] * Prob1[jj][a];
	      }
	      //accumulate score over individuals
	      if( h ){//heterozygous - prob H under null
		score[locus] -= 1.0 - H; 
	      }
	      else{//homozygous - prob (1-H) under null
		score[locus] += H;  // 0.5 * ( H  / ( 1.0 - H ) ); 
	      }
	      suminfo[locus] += H * (1.0 - H); // += 0.25 *( H /  (1.0 - H) );
	    }
	    ++locus;
	  }
	  //reset pointers ready to reuse next time	
	  free_matrix(Prob0, Loci->getNumberOfLoci(complocus));
	  free_matrix(Prob1, Loci->getNumberOfLoci(complocus));
	  // }
	++complocus;
      }
    }
    ind = 0;
  }
  for(int j = 0; j < NumLoci; ++j){
    sumscore[j] += score[j];
    sumscore2[j] += score[j] * score[j];
  }
  ++samples;
}

void HWTest::Output(const Vector_s LocusLabels){
  //header line
  outputfile <<"Locus\tScore\tCompleteInfo\tMissingInfo\tObservedInfo\tPercentInfo\tz-score\tp-value"<<endl;

 double EU, missing, complete, zscore;
  for(int j = 0; j < NumLoci; j++ ){
  //output locus labels from locus file
    //need same code as in ScoreTests to do for comp loci
    outputfile << LocusLabels[j] << "\t";

    EU = sumscore[ j ] / (double) samples;
    missing = sumscore2[ j ] / (double) samples - EU * EU;
    complete =  suminfo[ j ] / (double) samples;
    zscore = EU / sqrt( complete - missing );
    
    //outputfile.precision(2);
      outputfile << double2R(EU)                                << "\t"
		 << double2R(complete)                          << "\t"
		 << double2R(missing)                          << "\t"
		 << double2R(complete - missing)                << "\t"
		 << double2R(100*(complete - missing)/complete) << "\t"
		 << double2R(zscore)   << "\t"
		 << 2.0 * gsl_cdf_ugaussian_P (-fabs(zscore)) << "\t" << endl;//p-value
  }
}

string HWTest::double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << floor(x*100+0.5)/100.0;//for two decimal places
    return( ret.str() );
  }
}
