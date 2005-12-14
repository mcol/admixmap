/** 
 *   ADMIXMAP
 *   admixmap.cc 
 *   Top-level source file
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "admixmap.h"
#include <fstream>

using namespace std;

int ReadArgsFromFile(char* filename, int* xargc, char **xargv);
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, const string* const PopulationLabels);
void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      const Genome *Loci, Chromosome **Chrm, LogWriter& Log, const double *LogL, double coolness, bool anneal);
void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      LogWriter& Log);
void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(){

  cout << "\n-----------------------------------------------" << endl;
  cout << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  cout << "-----------------------------------------------" << endl;
  cout << "Programme Authors: " <<endl;
  cout << "David O'Donnell, Clive Hoggart and Paul McKeigue"<<endl;
  cout << "Copyright(c) 2002, 2003, 2004, 2005 LSHTM" <<endl;
  cout << "Send any comments or queries to david.odonnell @ucd.ie"<<endl;
  cout << "-----------------------------------------------"<<endl;
  cout << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl;
  cout << "under the terms of the GNU General Public License" <<endl;
  cout << "-----------------------------------------------" << endl;
}

int main( int argc , char** argv ){
  int    xargc = argc;
  char **xargv = argv;    

  if (argc < 2) {
    cout << "Please specify an options file or command-line arguments\n"
	 << "Usage:\n"
	 << "1. (not recommended) admixmap --[optionname]=[value] ...\n"
	 << "2. admixmap [optionfile], where optionfile is a text file containg a list of user options\n"
	 << "3. use a perl script. See sample perl script supplied with this program.\n"
	 << "Consult the manual for a list of user options."
	 << endl;
    exit(1); 
  } else if (argc == 2) {     // using options text file        
    xargc = 1;//NB initialise to 1 to mimic argc (arg 0 is prog name), otherwise first option is ignored later
    xargv = new char*[50];  // change 50 to max number of options
    ReadArgsFromFile(argv[1], &xargc, xargv);        
  }

  // ******************* PRIMARY INITIALIZATION ********************************************************************************
  //read user options
  AdmixOptions options(xargc, xargv);

  if(options.getDisplayLevel()>0)
    PrintCopyrightNotice();

  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  Log.StartMessage();

  smyrand( options.getSeed() );  // Initialise random number seed

  InputData data; //read data files and check (except allelefreq files)
  data.readData(&options, Log);//also sets 'numberofoutcomes' and 'populations' options

  //check user options
  options.checkOptions(Log, data.getNumberOfIndividuals());

  //print user options to args.txt; must be done after all options are set
  options.PrintOptions();

  Genome Loci;
  Loci.loadAlleleStatesAndDistances(&options, &data);//reads locusfile and creates CompositeLocus objects
  
  AlleleFreqs A(&Loci);
  A.Initialise(&options, &data, Log); //checks allelefreq files, initialises allele frequencies and finishes setting up Composite Loci
 
  Chromosome **chrm = 0; //Note: array of pointers to Chromosome
  chrm = Loci.GetChromosomes(options.getPopulations());  //create Chromosome objects
  Loci.SetSizes(Log);//prints length of genome, num loci, num chromosomes
    
  IndividualCollection *IC = new IndividualCollection(&options, &data, Loci, chrm);//NB call after A Initialise
  IC->LoadData(&options, &data);                             //and before L and R Initialise

  Latent L( &options, &Loci);    
  L.Initialise(IC->getSize(), data.GetPopLabels(), Log);

  Regression R[2];
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
    R[r].Initialise(r, IC, Log);
  Regression::OpenOutputFile(&options, IC, data.GetPopLabels(), Log);  

  if( options.isGlobalRho() )
    for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
      chrm[j]->InitialiseLociCorr(L.getrho());
    }
  IC->Initialise(&options, &Loci, data.GetPopLabels(), L.getrhoalpha(), L.getrhobeta(), Log, data.getMLEMatrix());
  //set expected Outcome
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
    R[r].SetExpectedY(IC);

  //  ******** single individual, one population, fixed allele frequencies  ***************************
  if( IC->getSize() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
    IC->getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
  // **************************************************************************************************
  else
    {
      int samples = options.getTotalSamples();
      // ******************* INITIALIZE TEST OBJECTS and ergodicaveragefile *******************************
      DispersionTest DispTest;
      StratificationTest StratTest;
      ScoreTests Scoretest;
      MisSpecAlleleFreqTest AlleleFreqTest;
      HWTest HWtest;
      std::ofstream avgstream; //output to ErgodicAverageFile
      
      if( options.getTestForDispersion() ){
	DispTest.Initialise(&options, Log, Loci.GetNumberOfCompositeLoci());    
      }
      if( options.getStratificationTest() )
	StratTest.Initialize( &options, Loci, chrm, IC, Log);
      if( options.getScoreTestIndicator() )
	Scoretest.Initialise(&options, IC, &Loci, chrm,data.GetPopLabels(), Log);
      if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	AlleleFreqTest.Initialise(&options, &Loci, Log );  
      if( options.getHWTestIndicator() )
	HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), Log);

      InitializeErgodicAvgFile(&options, IC, Log, &avgstream,data.GetPopLabels());

      //Write Initial values
      if(options.getIndAdmixHierIndicator()  ){
	if(options.getDisplayLevel()>2)Log.setDisplayMode(Quiet);
	else Log.setDisplayMode(Off);
	Log << "InitialParameterValues:\n";
	OutputParameters(-1, IC, &L, &A, R, &options, Log);
	Log << "\n";
      }

      string s = options.getResultsDir()+"/loglikelihoodfile.txt";
      ofstream loglikelihoodfile(s.c_str());

      // ******************* initialise stuff for annealing ************************************************
      bool anneal = (bool)(options.getAnnealIndicator()>0);
      double SumLogL = 0.0, SumLogLSq = 0.0, marg_L = 0.0;//modified log-likelihood, square and marginal likelihood
      double coolness = 1.0;
      std::ofstream annealstream;//for monitoring modified log likelihood when annealing
      int ci = options.getNumberOfAnnealedRuns();//defaults to starting at coolness of 1

      if(anneal){
	string s = options.getResultsDir()+"/annealmon.txt";
	annealstream.open(s.c_str());
	annealstream << "Coolness  Mean  Variance  logMargLikelihood"<<endl;
	ci = 0;
      }
      for( ;ci <= options.getNumberOfAnnealedRuns() ;ci++ ){//loop over temperature
	coolness = ci / (double)options.getNumberOfAnnealedRuns();
	//resets for start of each run
	SumLogL = 0.0;//cumulative sum of modified loglikelihood
	SumLogLSq = 0.0;//cumulative sum of square of modified loglikelihood
	if(anneal && options.getDisplayLevel() >0){
	  cout<<"\rSampling at coolness of "<<coolness;
	  if(!(ci%10) && ci < options.getNumberOfAnnealedRuns())cout<<"0";
	  if(ci == options.getNumberOfAnnealedRuns()) {
	    anneal = false;//finished annealed runs, last run is with unannealed likelihood
	    cout<<".00"<<endl;
	  }
	  cout<<flush;
	}
	if(options.getAnnealIndicator()==1)
	  Chromosome::setCoolness(coolness);//pass current coolness to Chromosome
      
	// *************************** BEGIN MAIN LOOP ******************************************************
	double LogL = IC->getLogLikelihood(&options, chrm, R, false);
	for( int iteration = 0; iteration <= samples; iteration++ ){
	  if(!anneal &&  !(iteration % options.getSampleEvery()) ){
	    WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
	  }

	  UpdateParameters(iteration, IC, &L, &A, R, &options, &Loci, chrm, Log, &LogL, coolness, anneal);

	  //compute loglikelihood and write to file
	  LogL = IC->getLogLikelihood(&options, chrm, R, (!anneal && iteration > options.getBurnIn()) );

	  if(!anneal)loglikelihoodfile<< iteration<<" " <<LogL<<endl;
	  //compute modified loglikelihood (at coolness of 1.0)
	  double lmod = LogL;
	  //if(coolness < 1.0){
	    lmod = IC->getModifiedLogLikelihood(&options, chrm, coolness);
	    //}
	  if(iteration > options.getBurnIn()){
	    SumLogL += lmod;
	    SumLogLSq += lmod;
	  }

	  if(!anneal){
	    // output every 'getSampleEvery()' iterations
	    if(!(iteration % options.getSampleEvery()) )
	      OutputParameters(iteration, IC, &L, &A, R, &options, Log);

	    // ** set merged haplotypes for allelic association score test 
	    if( iteration == options.getBurnIn() && options.getTestForAllelicAssociation() ){
	      Scoretest.SetAllelicAssociationTest(L.getalpha0());
	    }

	    //Updates and Output after BurnIn     
	    if( iteration > options.getBurnIn() ){
	      //dispersion test
	      if( options.getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(&A);
	      //stratification test
	      if( options.getStratificationTest() )StratTest.calculate(IC, A.GetAlleleFreqs(), Loci.GetChrmAndLocus(), 
								       options.getPopulations());
	      //score tests
	      if( options.getScoreTestIndicator() )
		Scoretest.Update(R[0].getDispersion());//possible error? what if 2 regression models?
	      //tests for mis-specified allelefreqs
	      if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
		AlleleFreqTest.Update(IC, &A, &Loci);
	      //test for Hardy-Weinberg eq
	      if( options.getHWTestIndicator() )
		HWtest.Update(IC, chrm, &Loci);
	  
	      // output every 'getSampleEvery() * 10' iterations (still after BurnIn)
	      if (!(iteration % (options.getSampleEvery() * 10))){    
	    
		//Ergodic averages
		if ( strlen( options.getErgodicAverageFilename() ) ){
		  int samples = iteration - options.getBurnIn();
		  if( options.getIndAdmixHierIndicator() ){
		    L.OutputErgodicAvg(samples,&avgstream);
		    for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
		      R[r].OutputErgodicAvg(samples, &avgstream);
		    A.OutputErgodicAvg(samples, &avgstream);
		  }
		  //if( IC->getSize()==1 )
		  IC->OutputErgodicAvg(samples, options.getMLIndicator(), &avgstream);
	      
		  avgstream << endl;
		}
		//Score Test output
		if( options.getScoreTestIndicator() )  Scoretest.Output(iteration,data.GetPopLabels());
	      }//end of 'every'*10 output
	    }//end if after BurnIn
	  }
	}//end main loop
	// *************************** END MAIN LOOP ******************************************************

	if(options.getAnnealIndicator()){//cannot use 'anneal' as it is false for final run
	  //calculate mean and variance of L_mod and write to file
	  double E_L_mod = 0.0,Var_L_mod = 0.0;
	  E_L_mod = SumLogL / ((double)options.getTotalSamples()-options.getBurnIn());
	  Var_L_mod = SumLogLSq/ ((double)options.getTotalSamples()-options.getBurnIn()) - E_L_mod * E_L_mod;
	  annealstream << coolness << "  "<<E_L_mod << "  " << Var_L_mod;
	  //use Simpson's rule to approximate integral
	  if( (ci == 0) || (ci == options.getNumberOfAnnealedRuns())) marg_L += E_L_mod;
	  else if( (ci%2) ) marg_L += 2.0*E_L_mod;
	  else marg_L += 4.0*E_L_mod; 
 
	  annealstream <<"  "<<marg_L / (double)(3*options.getNumberOfAnnealedRuns())<<endl;
	}
      }
      // *************************** END ANNEALING LOOP ******************************************************

      // *************************** OUTPUT AT END ***********************************************************
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
      else Log.setDisplayMode(On);
       if( options.getMLIndicator()){
	  IC->OutputChibEstimates(Log, options.getPopulations());
	if(IC->getSize()==1)
	  //MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
	  IC->OutputChibResults(Log);
	else{
	  ;
	}
      }
      IC->OutputDeviance(&options, chrm, R, Log, L.getSumLogRho(), Loci.GetNumberOfChromosomes());

      if(options.getAnnealIndicator()){
	Log<< "EstimatedLogMarginalLikelihood: "
	   <<  marg_L / (double)(3*options.getNumberOfAnnealedRuns())
	//<< "\nwith standard error of " 
	<< "\n";
      }

      //Residuals
      if(options.getNumberOfOutcomes() > 0)
	IC->OutputResiduals(options.getResidualFilename(), data.getOutcomeLabels(), options.getTotalSamples()-options.getBurnIn());
      //FST
      if( strlen( options.getHistoricalAlleleFreqFilename() ) ){
	A.OutputFST();
      }
      //stratification test
      if( options.getStratificationTest() ) StratTest.Output(Log);
      //dispersion test
      if( options.getTestForDispersion() )  DispTest.Output(options.getTotalSamples() - options.getBurnIn(), Loci, data.GetPopLabels());
      //tests for mis-specified allele frequencies
      if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	AlleleFreqTest.Output(options.getTotalSamples() - options.getBurnIn(), &Loci, data.GetPopLabels()); 
      //test for H-W eq
      if( options.getHWTestIndicator() )
	HWtest.Output(data.getLocusData()); 
      //finish writing score test output as R objects
      if( options.getScoreTestIndicator() ) Scoretest.ROutput();

      //output to likelihood ratio file
      if(options.getTestForAffectedsOnly())
	Individual::OutputLikRatios(options.getLikRatioFilename(), options.getTotalSamples()-options.getBurnIn(), data.GetPopLabels());
    
      if(annealstream.is_open())annealstream.close();
      if(avgstream.is_open())avgstream.close();
    }//end else

  // *************************** CLEAN UP ******************************************************  
  for(unsigned i = 0; i < Loci.GetNumberOfChromosomes(); i++){
    delete chrm[i];
  }
  A.CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), data.GetPopLabels());  
  delete IC;//must call explicitly so IndAdmixOutputter destructor finishes writing to indadmixture.txt
  delete []chrm;
   
  // ******************* acceptance rates - output to screen and log ***************************
  if( options.getIndAdmixHierIndicator() ){
    if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
    else Log.setDisplayMode(On);
#if POPADMIXSAMPLER == 2 
    if(options.getPopulations() > 1){
      Log << "Expected acceptance rate in admixture dispersion parameter sampler: "
	  << L.getEtaSamplerAcceptanceRate()
	  << "\nwith final step size of "
	  << L.getEtaSamplerStepsize() << "\n";
    }
#elif POPADMIXSAMPLER == 3
    if(options.getPopulations() > 1){
      Log << "Expected acceptance rate in admixture parameter Hamiltonian sampler: "
	  << L.getAlphaSamplerAcceptanceRate()
	  << "\nwith final step size of "
	  << L.getAlphaSamplerStepsize() << "\n";
    }
#endif
    
    if( options.isGlobalRho() && options.getPopulations() > 1 ){
      Log << "Expected acceptance rate in global sumintensities sampler: "
	  << L.getRhoSamplerAccRate()
	  << "\nwith final step size of "
	  << L.getRhoSamplerStepsize()
	  << "\n";
    }
    if(options.getCorrelatedAlleleFreqs()){
      Log<< "Expected acceptance rates in sampler for allele frequency proportion parameters: \n";
      for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	if(Loci(i)->GetNumberOfStates()>2)
	  Log << A.getAlphaSamplerAcceptanceRate(i) << " ";
      }
      Log<< "Expected acceptance rate in sampler for allele frequency dispersion parameter: \n";
      Log << A.getEtaRWSamplerAcceptanceRate(0)
	  << "\nwith final step sizes of \n";
      for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	if(Loci(i)->GetNumberOfStates()>2)
	  Log << A.getAlphaSamplerStepsize(i) << " " ;
      }
      Log <<  A.getEtaRWSamplerStepsize(0) << "\n" ;
    }
     
#if ETASAMPLER ==1
    if( strlen( options.getHistoricalAlleleFreqFilename() )){
      Log << "Expected acceptance rates in allele frequency dispersion parameter samplers:\n ";
      for(int k = 0; k < options.getPopulations(); ++k){Log << A.getEtaRWSamplerAcceptanceRate(k)<< " " ;}
      Log << "\nwith final step sizes of ";
      for(int k = 0; k < options.getPopulations(); ++k){Log <<  A.getEtaRWSamplerStepsize(k) << " ";}
      Log << "\n";
    }
#endif
  }
 
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  Log.ProcessingTime();
  if(options.getDisplayLevel()>0){
    cout << "Finished\n";
    for(unsigned i = 0; i < 80; ++i)cout<<"*";
    cout<<endl;
  }
  return 0;
}//end of main

int ReadArgsFromFile(char* filename, int* xargc, char **xargv){
  int  _maxLineLength=1024;
  ifstream fin(filename);
  std::string str;
  //read in line from file
  while (getline(fin,str,'\n')){

    if(str.length()){    //skip blank lines. **should skip also lines with only whitespace
      str.erase(0,str.find_first_not_of(" \t\n\r"));//trim leading whitespace
      if(str[0]!= '#'){ //ignore lines commented out with #
	if(str.find_first_of("#")<str.length())str.erase(str.find_first_of("#"));//ignore #comments
	//trim remaining whitespace
	str.erase(str.find_last_not_of(" \t\n\r")+1);//trailing whitespace
	if(str.find_first_of(" \t\n\r") <= str.length()){//check for any whitespace left
	  if(str.find_first_of(" \t\n\r") < str.find("="))//check for space before '='
	    str.erase(str.find_first_of(" \t\n\r"),str.find("=")-str.find_first_of(" \t\n\r"));//before '='
	  str.erase(str.find_first_of(" \t\n\r"),str.find_last_of(" \t\n")-str.find_first_of(" \t\n\r")+1);//after '='
	}
	//add line to xargv
	xargv[*xargc]=new char[_maxLineLength];
	strcpy(xargv[*xargc],"--");
	strcat(xargv[*xargc],str.c_str());
	++(*xargc);
      }}}
  fin.close();
  return 1;
}

//this function is here because three different objects have to write to avgstream
// for compatibility with parallelization, should rearrange output so that each object (L, R, A) writes one file 
// containing draws and ergodic averages for the parameters that it updates
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, const std::string* const PopulationLabels){
  Log.setDisplayMode(Quiet);
  //Open ErgodicAverageFile  
  if ( strlen( options->getErgodicAverageFilename() ) )
    {
      avgstream->open( options->getErgodicAverageFilename(), ios::out );
      if( !*avgstream ){
	Log.setDisplayMode(On);
	Log << "ERROR: Couldn't open Ergodic Average file\n";
	//exit( 1 );
      }
      else{
	Log << "Writing ergodic averages of parameters to "
	    << options->getErgodicAverageFilename() << "\n\n";
      }
      
      // Header line of ergodicaveragefile
      if( options->getIndAdmixHierIndicator() ){
	for( int i = 0; i < options->getPopulations(); i++ ){
	  *avgstream << "\""<<PopulationLabels[i] << "\" ";
	}
	if( options->isGlobalRho() )
	  *avgstream << " \"sumIntensities\"";
	else
	  *avgstream << "\"sumIntensities.mean\" ";
	
	
	// Regression parameters
	if( options->getNumberOfOutcomes() > 0 ){
	  for(int r = 0; r < individuals->getNumberOfOutcomeVars(); ++r){
	    *avgstream << "       \"intercept\" ";
	    if(strlen(options->getCovariatesFilename()) > 0){//if covariatesfile specified
	      for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
		*avgstream << individuals->getCovariateLabels(i) << " ";
	      }
	    }
	    if( !options->getTestForAdmixtureAssociation() ){
	      for( int k = 1; k < options->getPopulations(); k++ ){
		*avgstream << "\""<<PopulationLabels[k] << "\" ";
	      }
	    }
	    if( individuals->getOutcomeType(r)==0 )//linear regression
	      *avgstream << "       \"precision\" ";
	  }
	}
		
	// dispersion parameters
	if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	  for( int k = 0; k < options->getPopulations(); k++ ){
	    *avgstream << " \"eta" << k << "\" ";
	  }
	}
      }
      *avgstream << "\"MeanDeviance\"\t \"VarDeviance\"\t ";
      if(options->getMLIndicator()){//marginal likelihood calculation
	*avgstream<<"\"LogPrior\" \t \"LogPosterior\" \t \"LogPosteriorAdmixture\" \t \"LogPosteriorSumIntensities\" \t "
		  <<"\"LogPosteriorAlleleFreqs\" \t \"LogMarginalLikelihood\"";
      }
      *avgstream << "\n";
    }
  else
    {
      Log << "No ergodicaveragefile given\n";
    }
}

void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      const Genome *Loci, Chromosome **Chrm, LogWriter& Log, const double *LogL, double coolness, bool anneal){
  A->ResetAlleleCounts();
  
  // ** update global sumintensities
  if((options->getPopulations() > 1) && (IC->getSize() > 1) && 
     options->getIndAdmixHierIndicator() && (Loci->GetLengthOfGenome()> 0.0))
    L->UpdateRhoWithRW(IC, Chrm, *LogL);
  
  // ** Update individual-level parameters  
  IC->Update(iteration, options, Chrm, A, R, L->getpoptheta(), L->getalpha(), L->getrho(), L->getrhoalpha(), L->getrhobeta(),
	     Log, coolness, anneal);
  
  // ** update allele frequencies
  if(A->IsRandom()){
    A->Update((iteration > options->getBurnIn() && !anneal));
    IC->setGenotypeProbs(Chrm, Loci->GetNumberOfChromosomes());
    for(int i = 0; i < IC->getSize(); ++i)
      IC->getIndividual(i)->HMMIsBad(true); //if the allelefreqs are not fixed they are sampled between
    //individual updates. Therefore the forward probs in the HMMs must be updated and the current stored 
    //values of likelihood are invalid
  }
  //update population admixture Dirichlet parameters
  L->Update(iteration, IC, Log, anneal);
  
  // ** update regression parameters (if regression model)
  for(int r = 0; r < options->getNumberOfOutcomes(); ++r)
    R[r].Update((!anneal && iteration > options->getBurnIn()), IC);

}
void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      LogWriter& Log){
  if(options->getIndAdmixHierIndicator()  ){
    //Only output population-level parameters when there is a hierarchical model on indadmixture
    // ** pop admixture, sumintensities
    if(options->getPopulations() > 1) L->OutputParams(iteration,Log);
    //** dispersion parameter (if dispersion model)
	A->OutputEta(iteration, options, Log);
    // ** regression parameters
    for(int r = 0; r < options->getNumberOfOutcomes(); ++r)
      R[r].Output(iteration, options, Log);
    
    //** new line in logfile
	if( iteration == 0 ) {Log.setDisplayMode(Off);Log << "\n";}
  }
  if( options->getDisplayLevel()>2 ) cout << endl;
  if( iteration > options->getBurnIn() ){
    // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
    if ( strlen( options->getIndAdmixtureFilename() ) ) IC->OutputIndAdmixture();
    if(options->getOutputAlleleFreq())A->OutputAlleleFreqs();
  }
}
void WriteIterationNumber(const int iteration, const int width, int displayLevel){
  if( displayLevel > 2 ) {
    cout << setiosflags( ios::fixed );
    cout.width(width );
    cout << "\r"<< iteration << " ";
  }
  else if( displayLevel > 1 ){
    if(iteration==0)cout<<"Iterations so far:\n";
    cout << "\r"<< iteration<<flush;//displays iteration counter on screen
  }

}
