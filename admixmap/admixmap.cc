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
#include "Chromosome.h"//not needed if chrm is moved out
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

  /*----------------
  | Initialisation |
   ----------------*/
  //start timer
  long StartTime = time(0);
  tm timer;
  timer = *localtime( &StartTime );
  Log.StartMessage(&timer);

  options->checkOptions(&Log);
  // Initialise random number seed
  smyrand( options->getSeed() );

  //Initialise Objects
  InputData data;
  data.readData(options, &Log);

  Genome Loci;
  Loci.loadAlleleStatesAndDistances(options, &data, &Log);

  Chromosome **chrm = 0;
  IndividualCollection *IC = 0;
  chib MargLikelihood;

  std::vector<bool> _admixed;//don't belong
  bool _symmetric;          //here  
  Vector_d poptheta;
  std::string *PopulationLabels = 0;//possibly belongs in InputData

  AlleleFreqs A(&Loci);
  Latent L( options, &Loci, &Log);
  Regression R;

  A.LoadAlleleFreqs(options,&chrm,&Log,&data,&PopulationLabels);//NB this sets Populations option
  IC = new IndividualCollection(options,&data,Loci,chrm);//NB call after LoadAlleleFreqs
  IC->LoadData(options,&data, &Log);                             //and before L and R Initialise
 
  L.Initialise(IC, &LogFileStream, &_admixed, &_symmetric, &poptheta, PopulationLabels);
  if( options->getAnalysisTypeIndicator() >= 2)
    R.Initialise(IC, options, PopulationLabels, &Log);
  A.Initialise(options, data.getEtaPriorMatrix(), &Log, PopulationLabels);

  if( !options->getRhoIndicator() )
    for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
      chrm[j]->InitialiseLociCorr(L.getrho());
    }
  IC->Initialise(options, R.getbeta(), &Loci, PopulationLabels, L.getrhoalpha(), L.getrhobeta(), &Log, data.getMLEMatrix());

  options->PrintOptions();//NB: call after all options are set
                          //Currently all except Populations are set in AdmixOptions		
  /*-------------------------------------------------------
  |  single individual, one population, allele frequencies |
  ---------------------------------------------------------*/
  if( options->getAnalysisTypeIndicator() == -1 && options->getPopulations() == 1 && strlen(options->getAlleleFreqFilename()) )
    IC->getOnePopOneIndLogLikelihood(&Log, PopulationLabels, A.IsRandom());

  else{
    //initialise test objects
    DispersionTest DispTest;
    StratificationTest StratTest;
    ScoreTests Scoretest;
    MisSpecAlleleFreqTest AlleleFreqTest;
    HWTest HWtest;
 
    if( options->getTestForDispersion() ){
      DispTest.Initialise(options,&Log, A.GetNumberOfCompositeLoci());    
    }
    if( options->getStratificationTest() )
      StratTest.Initialize( options, Loci ,&Log);
    if( options->getScoreTestIndicator() )
      Scoretest.Initialise(options, IC, &Loci, chrm,PopulationLabels, &Log);
    if( options->getTestForMisspecifiedAlleleFreqs() || options->getTestForMisspecifiedAlleleFreqs2())
      AlleleFreqTest.Initialise(options, &Loci, &Log );  
    if( options->getHWTestIndicator() )
      HWtest.Initialise(options, Loci.GetTotalNumberOfLoci(), &Log);

    if( options->getTextIndicator() ){
      InitializeErgodicAvgFile(options,IC, &Log,&avgstream,PopulationLabels);
      }

 /*------------
  |  MAIN LOOP |
  ------------*/
    for( int iteration = 0; iteration <= options->getTotalSamples(); iteration++ ){
      if( !(iteration % options->getSampleEvery()) ){
	if( options->getAnalysisTypeIndicator() >= 0 && (!options->useCOUT() || iteration == 0) )
	  //output params to log when coutindicator = 0
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

      A.ResetAlleleCounts();
      //Update individual-level parameters  
      IC->Update(iteration, &A, &R, poptheta, options,
		 chrm, L.getalpha(), _symmetric, _admixed, L.getrhoalpha(), L.getrhobeta(),
		 &LogFileStream, &MargLikelihood);
      //update allele frequencies
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

      //update f summary for global rho
      if( !options->getRhoIndicator() )  
	for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
	  chrm[j]->SetLociCorr(L.getrho());
	}
      //update regression parameters (if regression model)
      if( options->getAnalysisTypeIndicator() >= 2)
	R.Update((iteration > options->getBurnIn()), IC);

      if( iteration == options->getBurnIn() && options->getTestForAllelicAssociation() ){
	A.SetMergedHaplotypes(L.getalpha0(), &LogFileStream, options->IsPedFile());
	Scoretest.SetAllelicAssociationTest();
      }
      
      // output every 'getSampleEvery()' iterations
      if( !(iteration % options->getSampleEvery()) ){
	if( options->getAnalysisTypeIndicator() >= 0 && options->getIndAdmixHierIndicator() ){
	  //Only output population-level parameters when there is a hierarchical model on indadmixture
	  //pop admixture, sumintensities
	  L.OutputParams(iteration, &LogFileStream);
	  //regression parameters
	  if( options->getAnalysisTypeIndicator() >= 2)
	    R.Output(iteration,&LogFileStream,options,IC);
	  //dispersion parameter (if dispersion model)
	  A.OutputEta(iteration, options, &LogFileStream);
	  //new line in logfile
	  if( !options->useCOUT() || iteration == 0 ) LogFileStream << endl;
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
	  Scoretest.Update(R.getDispersion(IC->getOutcomeType(0)));
	//tests for mis-specified allelefreqs
	if( options->getTestForMisspecifiedAlleleFreqs() || options->getTestForMisspecifiedAlleleFreqs2())
	  AlleleFreqTest.Update(IC, &A, &Loci);
	//test for Hardy-Weinberg eq
	if( options->getHWTestIndicator() )
	  HWtest.Update(IC, chrm, &Loci);

	// output every 'getSampleEvery() * 10' iterations (still after BurnIn)
	if (!(iteration % (options->getSampleEvery() * 10))){    
	  //FST
	  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	    A.OutputFST(options->IsPedFile());
	  }
	  //Ergodic averages
	  if ( strlen( options->getErgodicAverageFilename() ) ){
	    int samples = iteration - options->getBurnIn();
	    L.OutputErgodicAvg(samples,&avgstream);
	    if( options->getAnalysisTypeIndicator() >= 2)
	      R.OutputErgodicAvg(samples,IC,&avgstream);
	    A.OutputErgodicAvg(samples,options,&avgstream);
	    if( options->getAnalysisTypeIndicator() == -1 ){
	      IC->OutputErgodicAvg(samples,&MargLikelihood,&avgstream);
	    }
	    avgstream << endl;
	  }
	  //Test output

	  if( options->getStratificationTest() ) StratTest.Output();
	  if( options->getScoreTestIndicator() )
	    Scoretest.Output(iteration,PopulationLabels);
	}//end of 'every'*10 output
      }//end if after BurnIn


    }//end main loop

    //output at end
    //dispersion test
    if( options->getTestForDispersion() )  DispTest.Output(options->getTotalSamples() - options->getBurnIn(), Loci, PopulationLabels);
    //tests for mis-specified allele frequencies
    if( options->getTestForMisspecifiedAlleleFreqs() || options->getTestForMisspecifiedAlleleFreqs2())
      AlleleFreqTest.Output(options->getTotalSamples() - options->getBurnIn(), &Loci, PopulationLabels, options->IsPedFile()); 
    //test for H-W eq
   if( options->getHWTestIndicator() )
     HWtest.Output(options->IsPedFile(), data.getGeneInfoData()); 
    //Marginal Likelihood for a single individual
    if( options->getAnalysisTypeIndicator() == -1 )MargLikelihood.Output(&LogFileStream);
    //MLEs of admixture & sumintensities for nonhierarchical model on individual admixture
    if( options->getMLIndicator() )IC->Output(&LogFileStream);
    //finish writing score test output as R objects
    if( options->getScoreTestIndicator() ) Scoretest.ROutput();
  }//end else

  //  for(int i=0; i<A.getLoci()->GetNumberOfChromosomes(); i++){
  //     delete chrm[i];
  //   }

  A.CloseOutputFile(options->getTotalSamples() - options->getBurnIn(), PopulationLabels);  
  delete IC;//must call explicitly so IndAdmixOutputter destructor finishes writing to indadmixture.txt
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
    if( !options->getTestForAdmixtureAssociation() ){
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

