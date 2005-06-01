/*
  Class implements a score test for deviation from Hardy-Weinberg equilibrium
  in order to test for genotyping errors.
  This version evaluates for each single locus and sums over individuals.
  code for evaluation for each individual (rather than summing) is commented out.
*/
#include "HWTest.h"
#include "functions.h"

HWTest::HWTest(){
  sumscore = 0;
  sumscore2 = 0;
  suminfo = 0;
  NumInd = 0;
  NumLoci = 0;
}

HWTest::~HWTest(){
  //free_matrix(sumscore, NumInd);
  //free_matrix(sumscore2, NumInd);
  //free_matrix(suminfo, NumInd);

  delete[] sumscore;
  delete[] sumscore2;
  delete[] suminfo;
}

//void HWTest::Initialise(AdmixOptions *options, int nind, int nloci, LogWriter *Log){
void HWTest::Initialise(AdmixOptions *options, int nloci, LogWriter *Log){

  //NumInd = nind;
  NumLoci = nloci;

  if( options->getHWTestIndicator() ){//not really necessary
    if ( strlen( options->getHWTestFilename() ) ){
      outputfile.open( options->getHWTestFilename(), ios::out );
      if( !outputfile ){
	Log->logmsg(true,"ERROR: Couldn't open hwtestfile\n");
	exit( 1 );}
      else {
	Log->logmsg(true,"HW test file: ");    
	Log->logmsg(true,options->getHWTestFilename());
	Log->logmsg(true,"\n");
	
	//sumscore = alloc2D_d(NumInd, NumLoci);
	//sumscore2 = alloc2D_d(NumInd, NumLoci);
	//suminfo = alloc2D_d(NumInd, NumLoci);
	sumscore = new double[NumLoci];
	sumscore2 = new double[NumLoci];
	suminfo = new double[NumLoci];
      }
    }
    else{
      Log->logmsg(true,"No hwtestfile given\n");
      //exit(1);}
    }
  }
}    

void HWTest::Update(IndividualCollection *IC, Chromosome **C, Genome *Loci){
  double H, score;
  bool h;
  int locus = 0, complocus = 0;
  unsigned short **genotype = 0;
  Individual *ind = 0;
  double **Prob0 = 0, **Prob1 = 0;

  for(int i = 0; i < IC->getSize(); ++i){
    ind = IC->getIndividual(i);
    locus = 0;    //absolute locus number
    complocus = 0;//absolute comp locus number
    for(unsigned int chr = 0; chr < Loci->GetNumberOfChromosomes(); ++chr){ //loop over chromosomes
      for(unsigned int j = 0; j < C[chr]->GetSize(); ++j){                  //loop over comp loci on chromosome
	genotype = ind->getGenotype(complocus);

	//allocate arrays to hold marginal alleleprobs; could be done in function but easier to control here.      
	Prob0 = new double*[(*Loci)(complocus)->GetNumberOfLoci()];
      	Prob1 = new double*[(*Loci)(complocus)->GetNumberOfLoci()];

	for(int jj = 0; jj < (*Loci)(complocus)->GetNumberOfLoci(); ++jj){
	  Prob0[jj] = new double[ (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj)];
	  Prob1[jj] = new double[ (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj)];

	}
	//retrieve marginal alleleprobs for composite locus complocus, given current ancestry states on each gamete of ind
 	(*Loci)(complocus)->getLocusAlleleProbs(Prob0, ind->GetLocusAncestry( (int)chr, 0, j) );// #loci x #alleles array  
 	(*Loci)(complocus)->getLocusAlleleProbs(Prob1, ind->GetLocusAncestry( (int)chr, 1, j) );   

	H = 0.0;
	for(int jj = 0; jj < Loci->getNumberOfLoci(complocus); ++jj){       //loop over loci within comp locus
	  if( genotype[jj][0] != 0){ //non-missing genotype, assumes second gamete missing if first is
	    h = (genotype[jj][0] == genotype[jj][1]);

	    //compute prob of heterozygosity by summing component-wise products of off-diagonal marginal probs	    
	    for(int a0 = 0; a0 < (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj); ++a0)//loop over allele pairs
	      for(int a1 = 0; a1 < (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj); ++a1)
		H += (a0 != a1) * Prob0[jj][a0] * Prob1[jj][a1];
	    
	    score = 0.5 * ( -h + ( (1.0 - h) * H ) / ( 1.0 - H ) ); 
	    sumscore[locus] += score;
	    sumscore2[locus] += score * score;
	    suminfo[locus] += 0.25 * ( -h - ( ( 1.0 - h ) * H* H) / ( ( 1.0 - H ) * ( 1.0 - H ) ) );
	  }
	  ++locus;
	}
	//reset pointers ready to reuse next time	
	free_matrix(Prob0, Loci->getNumberOfLoci(complocus));
	free_matrix(Prob1, Loci->getNumberOfLoci(complocus));
	genotype = 0;
	++complocus;
      }
    }
    ind = 0;
  }
  ++samples;
}

void HWTest::Output(bool IsPedFile, Matrix_s LocusData){
 outputfile << "structure(.Data=c(" << endl;

  double EU, missing, complete;
  for(int j = 0; j < NumLoci; j++ ){
  //output locus labels from locus file
    //need same code as in ScoreTests to do for comp loci
      if(IsPedFile)
	outputfile << "\"" << LocusData[j+1][0] << "\"" << ",";
      else
	outputfile << LocusData[j][0] << ",";
      
      EU = sumscore[ j ] / (double) samples;
      missing = sumscore2[ j ] / (double) samples - EU * EU;
      complete =  suminfo[ j ] / (double) samples;
      
      outputfile << double2R(EU)                                << ",";
      outputfile << double2R(complete)                          << ",";
      outputfile << double2R(complete - missing)                << ",";
      outputfile << double2R(100*(complete - missing)/complete) << ",";
      outputfile << double2R(EU / sqrt( complete - missing ))   << "," << endl;
  }

    vector<int> dimensions(2,0);
    dimensions[0] = 6;
    dimensions[1] = NumLoci;

    vector<string> labels(6,"");
    labels[0] = "Locus";
    labels[1] = "Score";
    labels[2] = "CompleteInfo";
    labels[3] = "ObservedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "StdNormal";
    R_output3DarrayDimensions(&outputfile, dimensions, labels);
}

//copied from ScoreTests
void HWTest::R_output3DarrayDimensions(ofstream* stream,vector<int> dim,vector<string> labels)
{
  *stream << ")," << endl;
  *stream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    *stream << dim[i];
    if(i != dim.size() - 1){
      *stream << ",";
    }
  }
  *stream << ")," << endl;
  *stream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    *stream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      *stream << ",";
    }
  }
  *stream << "), character(0)))" << endl;
}

string HWTest::double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
