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
#include "IndividualCollection.h"
#include "Chromosome.h"
#include "chib.h"
#include "MisSpecAlleleFreqTest.h"
#include "HWTest.h"
#include <fstream>

using namespace std;

int ReadArgsFromFile(char* filename,int* xargc,char **xargv);
void InitializeErgodicAvgFile(AdmixOptions *options, IndividualCollection *individuals, LogWriter *Log, std::ofstream *avgstream, 
			      std::string *PopulationLabels);
void ProcessingTime(LogWriter*, long);

void PrintCopyrightNotice(){

    cout << "\n-----------------------------------------------" << endl;
    cout << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
    cout << "-----------------------------------------------" << endl;
    cout << "Programme Authors: " <<endl;
    cout << "David O'Donnell, Clive Hoggart, and Paul McKeigue"<<endl;
    cout << "Copyright(c) 2002, 2003, 2004, 2005 LSHTM" <<endl;
    cout << "Send any comments or queries to david.odonnell@ucd.ie"<<endl;
    cout << "-----------------------------------------------\n"<<endl;
    cout << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl;
    cout << "under the terms of the GNU General Public License" <<endl;
    cout << "-----------------------------------------------" << endl;
}

void submain(AdmixOptions* options){
  // This submain function is needed because the Latent object L is not destroyed properly if you use delete L. The ends of some R files which rely on a destructor being called are not generated without thus submain().
  //Not true any more - is submain still necessary?

  /*----------------
  | Initialisation |
   ----------------*/
  //start timer
  long StartTime = time(0);
  tm timer;
  timer = *localtime( &StartTime );

  LogWriter Log(options->getLogFilename(),options->useCOUT());
  Log.StartMessage(&timer);

  options->checkOptions(&Log);

  smyrand( options->getSeed() );  // Initialise random number seed

  InputData data; //read data files and check (except allelefreq files)
  data.readData(options, &Log);

  Genome Loci;
  Loci.loadAlleleStatesAndDistances(options, &data);//creates CompositeLocus objects
  
  AlleleFreqs A(&Loci);
  A.Initialise(options, &data, &Log); //checks allelefreq files, initialises allele frequencies and finishes setting up Composite Loci
  //Note: this sets Populations option

  options->PrintOptions();//NB: call after all options are set
                          //Currently all except Populations are set in AdmixOptions::SetOptions  
 
  Chromosome **chrm = 0; //Note: array of pointers to Chromosome
  chrm = Loci.GetChromosomes(options->getPopulations());  //create Chromosome objects
  Loci.SetSizes(&Log);//prints length of genome, num loci, num chromosomes
    
  IndividualCollection *IC = new IndividualCollection(options,&data,Loci,chrm);//NB call after A Initialise
  IC->LoadData(options,&data, &Log);                             //and before L and R Initialise

  Latent L( options, &Loci, &Log);    
  L.Initialise(IC->getSize(), data.GetPopLabels());

  Regression R0;
  Regression R1;
  if(options->getAnalysisTypeIndicator() > 1){
    R0.Initialise(0, IC, options, &Log);
    if(options->getAnalysisTypeIndicator() == 5)R1.Initialise(1, IC, options, &Log);
  }
  Regression::OpenOutputFile(options, IC, data.GetPopLabels(), &Log);  


  if( !options->getRhoIndicator() )
    for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
      chrm[j]->InitialiseLociCorr(L.getrho());
    }
  IC->Initialise(options, &Loci, data.GetPopLabels(), L.getrhoalpha(), L.getrhobeta(), &Log, data.getMLEMatrix());
  //set expected Y
  if(options->getAnalysisTypeIndicator() > 1){
    R0.SetExpectedY(IC);
    if(options->getAnalysisTypeIndicator() == 5)R1.SetExpectedY(IC);
  }		

  //   ** single individual, one population, allele frequencies 
   if( options->getAnalysisTypeIndicator() == -1 && options->getPopulations() == 1 && strlen(options->getAlleleFreqFilename()) )
     IC->getOnePopOneIndLogLikelihood(&Log, data.GetPopLabels());

  else
{
    //initialise test objects
    DispersionTest DispTest;
    StratificationTest StratTest;
    ScoreTests Scoretest;
    MisSpecAlleleFreqTest AlleleFreqTest;
    HWTest HWtest;
    std::ofstream avgstream; //output to ErgodicAverageFile
    chib MargLikelihood;
   
    if( options->getTestForDispersion() ){
      DispTest.Initialise(options,&Log, A.GetNumberOfCompositeLoci());    
    }
    if( options->getStratificationTest() )
      StratTest.Initialize( options, Loci, chrm, IC, &Log);
    if( options->getScoreTestIndicator() )
      Scoretest.Initialise(options, IC, &Loci, chrm,data.GetPopLabels(), &Log);
    if( options->getTestForMisspecifiedAlleleFreqs() || options->getTestForMisspecifiedAlleleFreqs2())
      AlleleFreqTest.Initialise(options, &Loci, &Log );  
    if( options->getHWTestIndicator() )
      HWtest.Initialise(options, Loci.GetTotalNumberOfLoci(), &Log);

    if( options->getTextIndicator() ){
      InitializeErgodicAvgFile(options,IC, &Log,&avgstream,data.GetPopLabels());
      }
    string s = options->getResultsDir()+"/loglikelihoodfile.txt";
    ofstream loglikelihoodfile(s.c_str());

 /*------------
  |  MAIN LOOP |
  ------------*/
    for( int iteration = 0; iteration <= options->getTotalSamples(); iteration++ ){
      if( !(iteration % options->getSampleEvery()) ){
	Log.Reset(iteration, (options->getAnalysisTypeIndicator() < 0), (int)( log10((double)options->getTotalSamples())+1 ) );
      }

      A.ResetAlleleCounts();

    //compute loglikelihood and write to file
    double LogL = 0.0;
    for(int i = 0; i < IC->getSize(); ++i)LogL += IC->getIndividual(i)->getLogLikelihood(options, chrm);
    loglikelihoodfile<< iteration<<" " <<LogL<<endl;

      // ** update global sumintensities
      if((options->getPopulations() > 1) && (IC->getSize() > 1) && options->getIndAdmixHierIndicator() && (Loci.GetLengthOfGenome()> 0.0))
	L.UpdateRhoWithRW(IC, chrm, LogL);
  
      // ** Update individual-level parameters  
      IC->Update(iteration, &A, &R0, &R1, L.getpoptheta(),options, chrm, L.getalpha(), L.getrhoalpha(), L.getrhobeta(),
		 &Log, &MargLikelihood);
      if((iteration %2)){
	L.Update(iteration, IC);//update pop admix params conditional on sums of ancestry states with jump indicators==1
	IC->ConjugateUpdateIndAdmixture(iteration, &R0, &R1, L.getpoptheta(),options, chrm, L.getalpha());//conjugate update of theta
      }

      // ** update allele frequencies
      A.Update((iteration > options->getBurnIn()));
      
      if( iteration > options->getBurnIn() ){
	if( options->getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(&A);
	if( options->getStratificationTest() )StratTest.calculate(IC, A.GetAlleleFreqs(), Loci.GetChrmAndLocus(), 
								  options->getPopulations());
      }  
 
      L.Update(iteration, IC);

      // ** update regression parameters (if regression model)
      if( options->getAnalysisTypeIndicator() >= 2){
	R0.Update((iteration > options->getBurnIn()), IC);
      if( options->getAnalysisTypeIndicator() == 5) R1.Update((iteration > options->getBurnIn()), IC);
      }

      // ** set merged haplotypes for allelic association score test 
     if( iteration == options->getBurnIn() && options->getTestForAllelicAssociation() ){
       Scoretest.SetAllelicAssociationTest(L.getalpha0());
      }
      
      // output every 'getSampleEvery()' iterations
      if( !(iteration % options->getSampleEvery()) ){
	if( options->getAnalysisTypeIndicator() >= 0 ){
	  //Only output population-level parameters when there is a hierarchical model on indadmixture
	  // ** pop admixture, sumintensities
	  L.OutputParams(iteration);
	  //** dispersion parameter (if dispersion model)
	  A.OutputEta(iteration, options, &Log);
	  // ** regression parameters
	  if( options->getAnalysisTypeIndicator() >= 2){
	    R0.Output(iteration, options, &Log);
	    if( options->getAnalysisTypeIndicator() == 5)
	      R1.Output(iteration, options, &Log);
	  }
	  //** new line in logfile
	  if( !options->useCOUT() || iteration == 0 ) Log.write("\n");
	}
	if( options->useCOUT() ) cout << endl;
	if( iteration > options->getBurnIn() ){
	  // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
	  if ( strlen( options->getIndAdmixtureFilename() ) ) IC->OutputIndAdmixture();
	  if(options->getOutputAlleleFreq())A.OutputAlleleFreqs();
	}
      }     
      //Updates and Output after BurnIn     
      if( iteration > options->getBurnIn() ){
	//score tests
	if( options->getScoreTestIndicator() )
	  Scoretest.Update(R0.getDispersion());//possible error? what if 2 regression models?
	//tests for mis-specified allelefreqs
	if( options->getTestForMisspecifiedAlleleFreqs() || options->getTestForMisspecifiedAlleleFreqs2())
	  AlleleFreqTest.Update(IC, &A, &Loci);
	//test for Hardy-Weinberg eq
	if( options->getHWTestIndicator() )
	  HWtest.Update(IC, chrm, &Loci);

	// output every 'getSampleEvery() * 10' iterations (still after BurnIn)
	if (!(iteration % (options->getSampleEvery() * 10))){    

	  //Ergodic averages
	  if ( strlen( options->getErgodicAverageFilename() ) ){
	    int samples = iteration - options->getBurnIn();
	    L.OutputErgodicAvg(samples,&avgstream);
	    if( options->getAnalysisTypeIndicator() >= 2){
	      R0.OutputErgodicAvg(samples, &avgstream);
	      if( options->getAnalysisTypeIndicator() == 5)
		R1.OutputErgodicAvg(samples, &avgstream);
	    }
	    A.OutputErgodicAvg(samples, &avgstream);
	    if( options->getAnalysisTypeIndicator() == -1 ){
	      IC->OutputErgodicAvg(samples,&MargLikelihood,&avgstream);
	    }
	    avgstream << endl;
	  }
	  //Score Test output
	  if( options->getScoreTestIndicator() )  Scoretest.Output(iteration,data.GetPopLabels());
	}//end of 'every'*10 output
      }//end if after BurnIn


    }//end main loop
    // ** output at end
    //FST
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      A.OutputFST();
    }
    //stratification test
    if( options->getStratificationTest() ) StratTest.Output(&Log);
    //dispersion test
    if( options->getTestForDispersion() )  DispTest.Output(options->getTotalSamples() - options->getBurnIn(), Loci, data.GetPopLabels());
    //tests for mis-specified allele frequencies
    if( options->getTestForMisspecifiedAlleleFreqs() || options->getTestForMisspecifiedAlleleFreqs2())
      AlleleFreqTest.Output(options->getTotalSamples() - options->getBurnIn(), &Loci, data.GetPopLabels()); 
    //test for H-W eq
    if( options->getHWTestIndicator() )
      HWtest.Output(data.getLocusData()); 
    //finish writing score test output as R objects
    if( options->getScoreTestIndicator() ) Scoretest.ROutput();

    //output to likelihood ratio file
    if(options->getTestForAffectedsOnly())
      Individual::OutputLikRatios(options->getLikRatioFilename(), options->getTotalSamples()-options->getBurnIn(), data.GetPopLabels());
    
    if( options->getMLIndicator() ){
      //MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
      IC->OutputChibEstimates(&Log, options->getPopulations());
      Log.logmsg(true,   "Log likelihood         (at estimates): ");Log.logmsg(true, MargLikelihood.getLogLikelihood());
      Log.logmsg(true, "\nLog prior              (at estimates): ");Log.logmsg(true, MargLikelihood.getLogPrior());
      Log.logmsg(true, "\nLog posterior          (at estimates): ");Log.logmsg(true, MargLikelihood.getLogPosterior());
      Log.logmsg(true, "\nLog marginal likelihood(at estimates): ");Log.logmsg(true, MargLikelihood.getLogMarginalLikelihood());
      Log.logmsg(true, "\n");
    }
    if(avgstream.is_open())avgstream.close();
 }//end else

  //  for(int i=0; i<A.getLoci()->GetNumberOfChromosomes(); i++){
  //     delete chrm[i];
  //   }
   A.CloseOutputFile((options->getTotalSamples() - options->getBurnIn())/options->getSampleEvery(), data.GetPopLabels());  
  delete IC;//must call explicitly so IndAdmixOutputter destructor finishes writing to indadmixture.txt
  delete []chrm;
  for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
    delete Loci(i);
  }

  // ** acceptance rates - output to screen and log
  if( options->getAnalysisTypeIndicator() >= 0 && options->getIndAdmixHierIndicator() ){
#if POPADMIXSAMPLER == 2 
    Log.logmsg(true,"Expected acceptance rate in admixture dispersion parameter sampler: ");
    Log.logmsg(true, L.getEtaSamplerAcceptanceRate());
    Log.logmsg(true, "\nwith final step size of ");
    Log.logmsg(true, L.getEtaSamplerStepsize());Log.logmsg(true, "\n");
    if(options->getPopulations() > 2){
      Log.logmsg(true,"Expected acceptance rate in admixture proportion parameter sampler: ");
      Log.logmsg(true, L.getMuSamplerAcceptanceRate());
      Log.logmsg(true, "\nwith final step size of ");
      Log.logmsg(true, L.getMuSamplerStepsize());Log.logmsg(true, "\n");
      }
#elif POPADMIXSAMPLER == 3 
    Log.logmsg(true,"Expected acceptance rate in admixture parameter Hamiltonian sampler: ");
    Log.logmsg(true, L.getAlphaSamplerAcceptanceRate());
    Log.logmsg(true, "\nwith final step size of ");
    Log.logmsg(true, L.getAlphaSamplerStepsize());Log.logmsg(true, "\n");
#endif
    
    if( !options->getRhoIndicator() ){
      Log.logmsg(true, "Expected acceptance rate in global sumintensities sampler: ");
      Log.logmsg(true, L.getRhoSamplerAccRate());
      Log.logmsg(true, "\nwith final step size of ");
      Log.logmsg(true, L.getRhoSamplerStepsize());
      Log.logmsg(true, "\n");
    }
    if(options->getCorrelatedAlleleFreqs()){
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
    if( strlen( options->getHistoricalAlleleFreqFilename() )){
      Log.logmsg(true, "Expected acceptance rates in allele frequency dispersion parameter samplers:\n ");
      for(int k = 0; k < options->getPopulations(); ++k){Log.logmsg(true, A.getEtaRWSamplerAcceptanceRate(k));Log.logmsg(true, " ");}
      Log.logmsg(true, "\nwith final step sizes of ");
      for(int k = 0; k < options->getPopulations(); ++k){Log.logmsg(true, A.getEtaRWSamplerStepsize(k));Log.logmsg(true, " ");}
      Log.logmsg(true, "\n");
    }
#endif
  }

  ProcessingTime(&Log, StartTime);
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

    AdmixOptions* options = new AdmixOptions;
    options->SetOptions(xargc, xargv);

    submain(options);
    cerr << "Finished\n";
  return 0;
}//end of main

int ReadArgsFromFile(char* filename,int* xargc,char **xargv){
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
void InitializeErgodicAvgFile(AdmixOptions *options, IndividualCollection *individuals, LogWriter *Log, std::ofstream *avgstream,
			      std::string *PopulationLabels){
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
    *avgstream << PopulationLabels[i] << " ";
  }
  if( !options->getRhoIndicator() )
    *avgstream << " \"sumIntensities\"";
  else
    *avgstream << "\"sumIntensities.beta\" ";


  // Regression parameters
  if( options->getAnalysisTypeIndicator() > 1 ){
    for(int r = 0; r < individuals->getNumberOfOutcomeVars(); ++r){
      *avgstream << "       \"intercept\" ";
      if(strlen(options->getCovariatesFilename()) > 0){//if covariatesfile specified
	for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
	  *avgstream << individuals->getCovariateLabels(i) << " ";
	}
      }
      if( !options->getTestForAdmixtureAssociation() ){
	for( int k = 1; k < options->getPopulations(); k++ ){
	  *avgstream << PopulationLabels[k] << " ";
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
  *avgstream << "\n";
    }
  else
    {
      Log->logmsg(true,"No ergodicaveragefile given\n");
    }
}

void ProcessingTime(LogWriter *Log, long StartTime)
{
  long Time = time(0);
  tm timer;
  timer = *localtime( &Time );
  

  Log->logmsg(false,"\nProgram finished at ");
  Log->logmsg(false,timer.tm_hour);
  Log->logmsg(false,":");
  Log->logmsg(false,timer.tm_min < 10 ? "0" : "");
  Log->logmsg(false,timer.tm_min);
  Log->logmsg(false,".");
  Log->logmsg(false,timer.tm_sec < 10 ? "0" : "" );
  Log->logmsg(false,timer.tm_sec);
  Log->logmsg(false," ");
  Log->logmsg(false,timer.tm_mday);
  Log->logmsg(false,"/");
  Log->logmsg(false,timer.tm_mon+1);
  Log->logmsg(false,"/");
  Log->logmsg(false,1900+timer.tm_year);
  Log->logmsg(false,"\n");

  Time -= StartTime;
  timer = *localtime(&Time);

  Log->logmsg(true,"Elapsed time = ");
  if( timer.tm_mday > 1 ){
    Log->logmsg(true,timer.tm_mday - 1 );
    Log->logmsg(true," day(s) ");
  }
    if( timer.tm_hour > 0 ){
  Log->logmsg(true,timer.tm_hour);
  Log->logmsg(true,"hour");if(timer.tm_hour > 1)Log->logmsg(true,"s");
    }
  Log->logmsg(true,timer.tm_min < 10 ? "0" : "");
  Log->logmsg(true,timer.tm_min);
  Log->logmsg(true,"m, ");
  Log->logmsg(true,timer.tm_sec < 10 ? "0" : "");
  Log->logmsg(true,timer.tm_sec);
  Log->logmsg(true,"s\n");

  //double realtime = difftime(time(0), StartTime);
//   realtime = pruntime();
//   Log->logmsg(true,"Elapsed time = ");
//   if(realtime > 3600.0){
//     Log->logmsg(true, (int)(realtime/3600));Log->logmsg(true,"hour");
//     if(realtime > 7200.0){Log->logmsg(true,"s");}
//     realtime = remainder(realtime, 3600.0);
//   }
//   if(realtime > 60.0){
//     Log->logmsg(true, (int)(realtime/3600));Log->logmsg(true,"min");
//     if(realtime > 120.0){Log->logmsg(true,"s");}
//     realtime = remainder(realtime, 60.0);
//   }
//   Log->logmsg(true, realtime);Log->logmsg(true, "seconds\n");
}

