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
			      LogWriter *Log, std::ofstream *avgstream, const string* const PopulationLabels);

void PrintCopyrightNotice(){

  cout << "\n-----------------------------------------------" << endl;
  cout << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  cout << "-----------------------------------------------" << endl;
  cout << "Programme Authors: " <<endl;
  cout << "David O'Donnell, Clive Hoggart and Paul McKeigue"<<endl;
  cout << "Copyright(c) 2002, 2003, 2004, 2005 LSHTM" <<endl;
  cout << "Send any comments or queries to david.odonnell@ucd.ie"<<endl;
  cout << "-----------------------------------------------"<<endl;
  cout << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl;
  cout << "under the terms of the GNU General Public License" <<endl;
  cout << "-----------------------------------------------" << endl;
}

int main( int argc , char** argv ){
  PrintCopyrightNotice();

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

  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(),options.useCOUT());
  Log.StartMessage();

  smyrand( options.getSeed() );  // Initialise random number seed

  InputData data; //read data files and check (except allelefreq files)
  data.readData(&options, &Log);//also sets 'numberofregressions' option

  //check user options
  options.checkOptions(&Log, data.getNumberOfIndividuals());

  Genome Loci;
  Loci.loadAlleleStatesAndDistances(&options, &data);//reads locusfile and creates CompositeLocus objects
  
  AlleleFreqs A(&Loci);
  A.Initialise(&options, &data, &Log); //checks allelefreq files, initialises allele frequencies and finishes setting up Composite Loci
  //Note: this sets Populations option

  //print user options to args.txt; must be done after all options are set
  options.PrintOptions();

  Chromosome **chrm = 0; //Note: array of pointers to Chromosome
  chrm = Loci.GetChromosomes(options.getPopulations());  //create Chromosome objects
  Loci.SetSizes(&Log);//prints length of genome, num loci, num chromosomes
    
  IndividualCollection *IC = new IndividualCollection(&options, &data, Loci, chrm);//NB call after A Initialise
  IC->LoadData(&options, &data);                             //and before L and R Initialise

  Latent L( &options, &Loci, &Log);    
  L.Initialise(IC->getSize(), data.GetPopLabels());

  Regression R[2];
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
    R[r].Initialise(r, IC, &Log);
  Regression::OpenOutputFile(&options, IC, data.GetPopLabels(), &Log);  

  if( options.isGlobalRho() )
    for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
      chrm[j]->InitialiseLociCorr(L.getrho());
    }
  IC->Initialise(&options, &Loci, data.GetPopLabels(), L.getrhoalpha(), L.getrhobeta(), &Log, data.getMLEMatrix());
  //set expected Outcome
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
    R[r].SetExpectedY(IC);

  //  ******** single individual, one population, fixed allele frequencies  ***************************
  if( IC->getSize() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
    IC->getOnePopOneIndLogLikelihood(&Log, data.GetPopLabels());
  // **************************************************************************************************
  else
    {

  // ******************* INITIALIZE TEST OBJECTS and ergodicaveragefile *******************************

      DispersionTest DispTest;
      StratificationTest StratTest;
      ScoreTests Scoretest;
      MisSpecAlleleFreqTest AlleleFreqTest;
      HWTest HWtest;
      std::ofstream avgstream; //output to ErgodicAverageFile
      
      if( options.getTestForDispersion() ){
	DispTest.Initialise(&options,&Log, Loci.GetNumberOfCompositeLoci());    
      }
      if( options.getStratificationTest() )
	StratTest.Initialize( &options, Loci, chrm, IC, &Log);
      if( options.getScoreTestIndicator() )
	Scoretest.Initialise(&options, IC, &Loci, chrm,data.GetPopLabels(), &Log);
      if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	AlleleFreqTest.Initialise(&options, &Loci, &Log );  
      if( options.getHWTestIndicator() )
	HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), &Log);

      InitializeErgodicAvgFile(&options,IC, &Log,&avgstream,data.GetPopLabels());

      string s = options.getResultsDir()+"/loglikelihoodfile.txt";
      ofstream loglikelihoodfile(s.c_str());

  // ******************* initialise stuff for annealing**************** *******************************
      bool anneal = options.getAnnealIndicator();
      double L_mod = 0.0, L_mod_sq = 0.0, marg_L = 0.0;//modified log-likelihood, square and marginal likelihood
      double coolness = 1.0;
      std::ofstream annealstream;//for monitoring modified log likelihood when annealing
      int ci = options.getNumberOfAnnealedRuns();//defaults to starting at coolness of 1

      if(anneal){
	string s = options.getResultsDir()+"/annealmon.txt";
	annealstream.open(s.c_str());
	annealstream << "Coolness  Mean  Variance  logMargLikelihood"<<endl;
	ci = 1;
      }
      for( ;ci <= options.getNumberOfAnnealedRuns() ;ci++ ){//loop over temperature
	coolness = (double)ci / (double)options.getNumberOfAnnealedRuns();
	//resets for start of each run
	if(anneal){
	  L_mod = 0.0;//ergodic sum of modified loglikelihood
	  L_mod_sq = 0.0;//ergodic sum of square of modified loglikelihood
	  cout<<"\rSampling at coolness of "<<coolness;
	  if(ci == options.getNumberOfAnnealedRuns()) {
	    anneal = false;//finished annealed runs, last run is with unannealed likelihood
	    cout<<".00"<<endl;
	  }
	}
	Chromosome::setCoolness(coolness, &L_mod);//pass current coolness and pointer to L_mod to Chromosome
      
	// *************************** BEGIN MAIN LOOP ******************************************************
	for( int iteration = 0; iteration <= options.getTotalSamples(); iteration++ ){
	  if(!anneal &&  !(iteration % options.getSampleEvery()) ){
	    if(IC->getSize() >1)
	      Log.Reset(iteration, (int)( log10((double)options.getTotalSamples())+1 ) );
	    if( options.useCOUT() ) {
	      cout << setiosflags( ios::fixed );
	      cout.width((int)( log10((double)options.getTotalSamples())+1 ) );
	      cout << iteration << " ";
	    }
	    else cout << "\r"<< iteration<<flush;//displays iteration counter on screen
	  }
	
	  A.ResetAlleleCounts();
	
	  //compute loglikelihood and write to file
	  double LogL = IC->getLogLikelihood(&options, chrm, options.getAnnealIndicator());
	  //includes calculation of unannealed loglikelihood
	  L_mod_sq += L_mod * L_mod;

	  if(!anneal)loglikelihoodfile<< iteration<<" " <<LogL<<endl;
	
	  // ** update global sumintensities
	  if((options.getPopulations() > 1) && (IC->getSize() > 1) && 
	     options.getIndAdmixHierIndicator() && (Loci.GetLengthOfGenome()> 0.0))
	    L.UpdateRhoWithRW(IC, chrm, LogL);
	
	  // ** Update individual-level parameters  
	  IC->Update(iteration, &options, chrm, &A, R, L.getpoptheta(), L.getalpha(), L.getrho(), L.getrhoalpha(), L.getrhobeta(),
		     &Log, anneal);
	  //if((iteration %2)){
	  //L.Update(iteration, IC);//update pop admix params conditional on sums of ancestry states with jump indicators==1
	  //IC->ConjugateUpdateIndAdmixture(iteration, R, L.getpoptheta(),options, chrm, L.getalpha());//conjugate update of theta
	  //       }
	
	  // ** update allele frequencies
	  if(A.IsRandom()){
	    A.Update((iteration > options.getBurnIn()));
	    for(int i = 0; i < IC->getSize(); ++i)
	      IC->getIndividual(i)->HMMIsBad(true); //if the allelefreqs are not fixed they are sampled between
	    //individual updates. Therefore the forward probs in the HMMs must be updated and the current stored 
	    //values of likelihood are invalid
	  }
	
	  if( !anneal && iteration > options.getBurnIn() ){
	    if( options.getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(&A);
	    if( options.getStratificationTest() )StratTest.calculate(IC, A.GetAlleleFreqs(), Loci.GetChrmAndLocus(), 
								      options.getPopulations());
	  }  
	
	  L.Update(iteration, IC);
	
	  // ** update regression parameters (if regression model)
	  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
	    R[r].Update((iteration > options.getBurnIn()), IC);
	
	  // ** set merged haplotypes for allelic association score test 
	  if( !anneal && iteration == options.getBurnIn() && options.getTestForAllelicAssociation() ){
	    Scoretest.SetAllelicAssociationTest(L.getalpha0());
	  }
	
	  // output every 'getSampleEvery()' iterations
	  if(!anneal &&  !(iteration % options.getSampleEvery()) ){
	    if( options.getIndAdmixHierIndicator() ){
	      //Only output population-level parameters when there is a hierarchical model on indadmixture
	      // ** pop admixture, sumintensities
	      if(options.getPopulations() > 1) L.OutputParams(iteration);
	      //** dispersion parameter (if dispersion model)
		  A.OutputEta(iteration, &options, &Log);
	      // ** regression parameters
	      for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
		R[r].Output(iteration, &options, &Log);
	    
	      //** new line in logfile
		  if( !options.useCOUT() || iteration == 0 ) Log.write("\n");
	    }
	    if( options.useCOUT() ) cout << endl;
	    if( iteration > options.getBurnIn() ){
	      // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
	      if ( strlen( options.getIndAdmixtureFilename() ) ) IC->OutputIndAdmixture();
	      if(options.getOutputAlleleFreq())A.OutputAlleleFreqs();
	    }
	  }     
	  //Updates and Output after BurnIn     
	  if( !anneal && iteration > options.getBurnIn() ){
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
		L.OutputErgodicAvg(samples,&avgstream);
		for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
		  R[r].OutputErgodicAvg(samples, &avgstream);
		A.OutputErgodicAvg(samples, &avgstream);
		//if( IC->getSize()==1 )
		IC->OutputErgodicAvg(samples, options.getMLIndicator(), &avgstream);
	      
		avgstream << endl;
	      }
	      //Score Test output
	      if( options.getScoreTestIndicator() )  Scoretest.Output(iteration,data.GetPopLabels());
	    }//end of 'every'*10 output
	  }//end if after BurnIn
	
	}//end main loop
	// *************************** END MAIN LOOP ******************************************************

	if(options.getAnnealIndicator()){//cannot use 'anneal' as it is false for final run
	  //calculate mean and variance of L_mod and write to file
	  double E_L_mod = 0.0,Var_L_mod = 0.0;
	  E_L_mod = L_mod / ((double)options.getTotalSamples()-options.getBurnIn() + 1.0);
	  Var_L_mod = L_mod_sq/ ((double)options.getTotalSamples()-options.getBurnIn() + 1.0) - E_L_mod * E_L_mod;
	  annealstream << coolness << "  "<<E_L_mod << "  " << Var_L_mod; 
	  marg_L += 0.01 * L_mod;
	  annealstream <<"  "<<marg_L / ((double)options.getTotalSamples()-options.getBurnIn() + 1.0)<<endl;
	}
      }
      // *************************** END ANNEALING LOOP ******************************************************

      // *************************** OUTPUT AT END ***********************************************************
      if( options.getMLIndicator() ){
	//MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
	IC->OutputChibEstimates(&Log, options.getPopulations());
      }
      IC->OutputDeviance(&options, chrm, R, &Log, L.getSumLogRho(), Loci.GetNumberOfChromosomes());

      if(options.getAnnealIndicator()){
	Log.logmsg(true, "Log Marginal Likelihood from simulated annealing: ");
	Log.logmsg(true, marg_L / ((double)options.getTotalSamples()-options.getBurnIn() + 1.0));
	//Log.logmsg(true, "\nwith standard error of ");Log.logmsg(true, );Log.logmsg(true, "\n");
      }

      //FST
      if( strlen( options.getHistoricalAlleleFreqFilename() ) ){
	A.OutputFST();
      }
      //stratification test
      if( options.getStratificationTest() ) StratTest.Output(&Log);
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
#if POPADMIXSAMPLER == 2 
    if(options.getPopulations() > 1){
      Log.logmsg(true,"Expected acceptance rate in admixture dispersion parameter sampler: ");
      Log.logmsg(true, L.getEtaSamplerAcceptanceRate());
      Log.logmsg(true, "\nwith final step size of ");
      Log.logmsg(true, L.getEtaSamplerStepsize());Log.logmsg(true, "\n");
      //     if(options.getPopulations() > 2){
      //       Log.logmsg(true,"Expected acceptance rate in admixture proportion parameter sampler: ");
      //       Log.logmsg(true, L.getMuSamplerAcceptanceRate());
      //       Log.logmsg(true, "\nwith final step size of ");
      //       Log.logmsg(true, L.getMuSamplerStepsize());Log.logmsg(true, "\n");
      //       }
    }
#elif POPADMIXSAMPLER == 3
    if(options.getPopulations() > 1){
      Log.logmsg(true,"Expected acceptance rate in admixture parameter Hamiltonian sampler: ");
      Log.logmsg(true, L.getAlphaSamplerAcceptanceRate());
      Log.logmsg(true, "\nwith final step size of ");
      Log.logmsg(true, L.getAlphaSamplerStepsize());Log.logmsg(true, "\n");
    }
#endif
    
    if( options.isGlobalRho() && options.getPopulations() > 1 ){
      Log.logmsg(true, "Expected acceptance rate in global sumintensities sampler: ");
      Log.logmsg(true, L.getRhoSamplerAccRate());
      Log.logmsg(true, "\nwith final step size of ");
      Log.logmsg(true, L.getRhoSamplerStepsize());
      Log.logmsg(true, "\n");
    }
    if(options.getCorrelatedAlleleFreqs()){
      Log.logmsg(true, "Expected acceptance rates in sampler for allele frequency prior parameters: \n");
      for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	if(Loci(i)->GetNumberOfStates()>2)
	  Log.logmsg(true, A.getAlphaSamplerAcceptanceRate(i));Log.logmsg(true, " ");
      }
      Log.logmsg(true, A.getEtaRWSamplerAcceptanceRate(0));
      Log.logmsg(true, "\nwith final step sizes of \n");
      for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	if(Loci(i)->GetNumberOfStates()>2)
	  Log.logmsg(true, A.getAlphaSamplerStepsize(i));Log.logmsg(true, " ");
      }
      Log.logmsg(true, A.getEtaRWSamplerStepsize(0));
      Log.logmsg(true, "\n");
    }
     
#if ETASAMPLER ==1
    if( strlen( options.getHistoricalAlleleFreqFilename() )){
      Log.logmsg(true, "Expected acceptance rates in allele frequency dispersion parameter samplers:\n ");
      for(int k = 0; k < options.getPopulations(); ++k){Log.logmsg(true, A.getEtaRWSamplerAcceptanceRate(k));Log.logmsg(true, " ");}
      Log.logmsg(true, "\nwith final step sizes of ");
      for(int k = 0; k < options.getPopulations(); ++k){Log.logmsg(true, A.getEtaRWSamplerStepsize(k));Log.logmsg(true, " ");}
      Log.logmsg(true, "\n");
    }
#endif
  }
 
  Log.ProcessingTime();
  cout << "Finished\n";
  for(unsigned i = 0; i < 80; ++i)cout<<"*";
  cout<<endl;
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
			      LogWriter *Log, std::ofstream *avgstream, const std::string* const PopulationLabels){
  //Open ErgodicAverageFile  
  if ( strlen( options->getErgodicAverageFilename() ) )
    {
      avgstream->open( options->getErgodicAverageFilename(), ios::out );
      if( !*avgstream ){
	Log->logmsg(true,"ERROR: Couldn't open Ergodic Average file\n");
	//exit( 1 );
      }
      else{
	Log->logmsg(true,"Writing ergodic averages of parameters to ");
	Log->logmsg(true,options->getErgodicAverageFilename());
	Log->logmsg(true,"\n\n");
      }
      
      // Header line of ergodicaveragefile
      for( int i = 0; i < options->getPopulations(); i++ ){
	*avgstream << "\""<<PopulationLabels[i] << "\" ";
      }
      if( options->isGlobalRho() )
	*avgstream << " \"sumIntensities\"";
      else
	*avgstream << "\"sumIntensities.beta\" ";
      
      
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
	    *avgstream << "       \"precision\"";
	}
      }
      
      
      // dispersion parameters
      if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	for( int k = 0; k < options->getPopulations(); k++ ){
	  *avgstream << " \"eta" << k << "\"";
	}
      }
      *avgstream << "\"MeanDeviance\"\t \"VarDeviance\"\t ";
      if(options->getMLIndicator()){//marginal likelihood calculation
	*avgstream<<"\"LogMarginalLikelihood \" ";
      }
      *avgstream << "\n";
    }
  else
    {
      Log->logmsg(true,"No ergodicaveragefile given\n");
    }
}

