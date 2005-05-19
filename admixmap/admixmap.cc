
#include "admixmap.h"
#include "IndividualCollection.h"
#include "Chromosome.h"//not needed if chrm is moved out
#include "chib.h"
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
    cout << "Clive Hoggart, Richard Sharp, Nigel Wetters, David O'Donnell and Paul McKeigue"<<endl;
    cout << "Copyright(c) 2002, 2003, 2004, 2005 LSHTM" <<endl;
    cout << "Send any comments or queries to david.odonnell@ucd.ie"<<endl;
    cout << "-----------------------------------------------\n"<<endl;
    cout << "This program is free software; you can redistribute it and/or modify"<<endl;
    cout << "it under the terms of the GNU General Public License as published by"<<endl;
    cout << "the Free Software Foundation; either version 2 of the License, or"<<endl;
    cout << "any later version."<<endl;
    cout << "\n";
    cout << "This program is distributed in the hope that it will be useful,"<<endl;
    cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of"<<endl;
    cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"<<endl;
    cout << "GNU General Public License for more details."<<endl;
    cout << "\n";
    cout << "You should have received a copy of the GNU General Public License"<<endl;
    cout << "along with this program; if not, write to the Free Software"<<endl;
    cout << "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307"<<endl;
    cout << "USA or contact the email address above."<<endl;
    cout << "-----------------------------------------------" << endl;
}

void submain(AdmixOptions* options){
  // This submain function is needed because the Latent object L is not destroyed properly if you use delete L. The ends of some R files which rely on a destructor being called are not generated without thus submain().
  //Not true any more - is submain still necessary?

  std::ofstream LogFileStream;//output to logfile
  std::ofstream avgstream; //output to ErgodicAverageFile

  LogFileStream.open( options->getLogFilename(), ios::out ); 
  LogWriter Log;
  Log.Initialise(&LogFileStream,options->useCOUT());

  options->checkOptions(&Log);//should be in AdmixOptions  

  InputData data;
  data.readData(options, &Log);

  Genome Loci;
  Chromosome **chrm = 0;
  IndividualCollection *IC;
  IC = 0;
  StratificationTest StratTest;
  ScoreTests Scoretest;
  DispersionTest DispTest;
  chib MargLikelihood;

  std::vector<bool> _admixed;
  bool _symmetric;


  Vector_d poptheta;
  std::string *PopulationLabels = 0;//possibly belongs in InputData

  /*----------------
  | Initialisation |
   ----------------*/
  long StartTime = time(0);
  
  tm timer;
  timer = *localtime( &StartTime );
  Log.StartMessage(options->getTotalSamples(), options->getBurnIn(),&timer);

  // Initialise random number seed
  smyrand( options->getSeed() );
  //Initialise Objects
  AlleleFreqs A(&Loci);
  Latent L( options, &Loci, &Log);
  Regression R;

  A.LoadAlleleFreqs(options,&chrm,&Log,&data,&PopulationLabels);//NB this sets Populations option
  data.determineIfPedFile(options);
  IC = new IndividualCollection(options,&data,Loci,chrm);//NB call after LoadAlleleFreqs
  IC->LoadGenotypes(options,&data, &Log, &Loci);                             //and before L and R Initialise
 
  L.Initialise(IC, &LogFileStream, &_admixed, &_symmetric, &poptheta, PopulationLabels);
  R.Initialise(IC, options, PopulationLabels, &Log);
  A.Initialise(options, data.getEtaPriorMatrix(), &Log, PopulationLabels, L.getrho());
  IC->Initialise(options, R.getbeta(), &Loci, PopulationLabels, L.getrhoalpha(), L.getrhobeta(), &Log, data.getMLEMatrix());

  options->PrintOptions();//NB: call after all options are set
                          //Currently all except Populations are set in AdmixOptions		
  /*-------------------------------------------------------
  |  single individual, one population, allele frequencies |
  ---------------------------------------------------------*/
  if( options->getAnalysisTypeIndicator() == -1 && options->getPopulations() == 1 && strlen(options->getAlleleFreqFilename()) )
    IC->getOnePopOneIndLogLikelihood(&Log,&A,PopulationLabels);

  else{
   
    Scoretest.Initialise(options, IC, &Loci, chrm,PopulationLabels, &Log);  
    StratTest.Initialize( options, Loci ,&Log);
    DispTest.Initialise(options,&Log, A.GetNumberOfCompositeLoci());
    if( options->getTextIndicator() ){
      InitializeErgodicAvgFile(options,IC, &Log,&avgstream,PopulationLabels);
      }

 /*------------
  |  MAIN LOOP |
  ------------*/
    for( int iteration = 0; iteration <= options->getTotalSamples(); iteration++ ){
      if( !(iteration % options->getSampleEvery()) ){
	if( options->getAnalysisTypeIndicator() >= 0 && (!options->useCOUT() || iteration == 0) )
	  //do we really want to output pars to log when coutindicator = 0?
	  {
	    LogFileStream << setiosflags( ios::fixed );
	    LogFileStream.width( (int)( log10((double)options->getTotalSamples())+1 ) );
 	    LogFileStream << iteration << " ";
	  }
	if( options->useCOUT() ) {
	  cout << setiosflags( ios::fixed );
	  cout.width( (int)( log10((double)options->getTotalSamples())+1 ) );
	  cout << iteration << " ";
	}
      }
//Resets before updates
      A.Reset();
      //Updates  
      IC->Update(iteration, &A, &R, poptheta, options,
		 chrm, L.getalpha(), _symmetric, _admixed, L.getrhoalpha(), L.getrhobeta(),
		 &LogFileStream, &MargLikelihood);
      A.Update(iteration,options->getBurnIn());
      
      if( iteration > options->getBurnIn() ){
	if( options->getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(&A);
	if( options->getStratificationTest() )StratTest.calculate(IC, A.GetAlleleFreqs(), Loci.GetChrmAndLocus());
      }  
 
      // Latent should not need to know anything about the number or positions of loci
      // with a global rho model, update of rho should be via a Metropolis random walk conditioned on the HMM likelihood
      // with a hierarchical rho model, update of hyperparameters should be via sufficient statistics: 
      // sum of rho and rho-squared over all individuals or gametes 
      L.Update(iteration, IC, &poptheta,&LogFileStream);
      A.ResetSumAlleleFreqs();
      if( !options->getRhoIndicator() )  A.load_f(L.getrho(),chrm);
      R.Update(IC);
      
      if( iteration == options->getBurnIn() && options->getTestForAllelicAssociation() ){
	A.SetMergedHaplotypes(L.getalpha0(), &LogFileStream, options->IsPedFile());
	Scoretest.SetAllelicAssociationTest();
      }
      
      // output every 'getSampleEvery()' iterations
      if( !(iteration % options->getSampleEvery()) ){
	if( options->getAnalysisTypeIndicator() >= 0 && options->getIndAdmixHierIndicator() ){
	  //Only output population-level parameters when there is a hierarchical model on indadmixture
	  L.OutputParams(iteration, &LogFileStream);
	  R.Output(iteration,&LogFileStream,options,IC);
	  A.OutputEta(iteration, options, &LogFileStream);
	  if( !options->useCOUT() || iteration == 0 ) LogFileStream << endl;
	}
	if( options->useCOUT() ) cout << endl;
	if( iteration > options->getBurnIn() ){
	  // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
	  if ( strlen( options->getIndAdmixtureFilename() ) ) IC->OutputIndAdmixture();
	  if(options->getOutputAlleleFreq())A.OutputAlleleFreqs();
	}
      }     
      //Output and scoretest updates after BurnIn     
      if( iteration > options->getBurnIn() ){
	Scoretest.Update(R.getDispersion(IC->getOutcomeType(0)), &A);
	R.SumParameters(options->getAnalysisTypeIndicator());
	
	// output every 'getSampleEvery() * 10' iterations
	if (!(iteration % (options->getSampleEvery() * 10))){    
	  //FST
	  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	    A.OutputFST(options->IsPedFile());
	  }
	  //Ergodic averages
	  if ( strlen( options->getErgodicAverageFilename() ) ){
	    int samples = iteration - options->getBurnIn();
	    L.OutputErgodicAvg(samples,&avgstream);
	    R.OutputErgodicAvg(samples,IC,&avgstream);
	    A.OutputErgodicAvg(samples,options,&avgstream);
	    if( options->getAnalysisTypeIndicator() == -1 ){
	      IC->OutputErgodicAvg(samples,&MargLikelihood,&avgstream);
	    }
	    avgstream << endl;
	  }
	  //Test output
	  if( options->getTestForDispersion() )  DispTest.Output(iteration - options->getBurnIn(), Loci);
	  if( options->getStratificationTest() ) StratTest.Output();
	  Scoretest.Output(iteration,PopulationLabels);
	}//end of 'every'*10 output
      }//end if after BurnIn
    }//end main loop
    
    if( options->getAnalysisTypeIndicator() == -1 )MargLikelihood.Output(&LogFileStream);
    if( options->getMLIndicator() )IC->Output(&LogFileStream);
    Scoretest.ROutput();
  }//end else

  //  for(int i=0; i<A.getLoci()->GetNumberOfChromosomes(); i++){
  //     delete chrm[i];
  //   }
  
  delete IC;//must call explicitly as IndAdmixOutputter destructor finishes writing to indadmixture.txt
  delete []chrm;
  for(unsigned int i=0; i < Loci.GetNumberOfCompositeLoci(); i++){
    delete Loci(i);
  }
  //delete []PopulationLabels;
  
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
      Log->logmsg(true,"Writing ergodic averages to ");
      Log->logmsg(true,options->getErgodicAverageFilename());
      Log->logmsg(true,"\n");
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
    *avgstream << "       \"intercept\" ";
    if( individuals->GetNumberOfInputRows() == individuals->getSize() ){
      for( int i = 0; i < individuals->GetNumberOfInputCols(); i++ ){
	*avgstream << individuals->getCovariateLabels(i) << " ";
      }
    }
    if( !options->getScoreTestIndicator() ){
      for( int k = 1; k < options->getPopulations(); k++ ){
	*avgstream << PopulationLabels[k] << " ";
      }
    }
  }
  if( options->getAnalysisTypeIndicator() == 2 )
    *avgstream << "       \"precision\"";

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
  

  Log->logmsg(false,"Program finished at ");
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

