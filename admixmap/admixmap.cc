
#include "admixmap.h"
#include "IndividualCollection.h"
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
    cout << "Copyright(c) 2002, 2003, 2004 LSHTM" <<endl;
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
  // This submain function is needed because the Latent object L does not destruct properly if you use delete L. The ends of some R files which rely on a destructor being called are not generated without thus submain().

  std::ofstream LogFileStream;//output to logfile
  std::ofstream avgstream; //output to ErgodicAverageFile

  LogFileStream.open( options->getLogFilename(), ios::out ); 
  LogWriter Log;
  Log.Initialise(&LogFileStream,options->useCOUT());

  options->checkOptions(&Log);//should be in AdmixOptions  

  InputData data;
  data.readData(options, &Log);

  Genome *chrm;//doesn't belong here
  IndividualCollection *IC;
  IC = 0;
  StratificationTest StratTest;
  ScoreTests Scoretest;
  DispersionTest DispTest;
  chib MargLikelihood;

  std::vector<bool> _admixed;
  bool _symmetric;
  // LociCorrSummary is a vector of terms of the form exp( - rho*x_i) where x_i is the map distance between two adjacent loci
  // with a global rho model, this vector is same for all individuals and calculated only once. 
  // should be calculated and stored in the Loci object
  Vector_d LociCorrSummary;//summary of correlation in ancestry between loci,was called f in Latent
  Vector_d SumLogTheta;
  Vector_d poptheta;
  std::string *PopulationLabels;//possibly belongs in InputData

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
  AlleleFreqs A;
  Latent L( options, A.getLoci(), &Log);
  Regression R;

  A.LoadAlleleFreqs(options,&chrm,&Log,&data,&PopulationLabels);//NB this sets Populations option
  IC = new IndividualCollection(options,data.getGeneticData(),*(A.getLoci()),*chrm);//NB call after LoadAlleleFreqs
  IC->LoadGenotypes(options,&data, &Log, A.getLoci());                             //and before L and R Initialise

  L.Initialise(IC,&LogFileStream, &_admixed,&_symmetric,&LociCorrSummary,&poptheta);
  R.Initialise(IC,options, &Log);
  A.Initialise(options,data.getEtaPriorMatrix(),&Log,PopulationLabels);
  IC->Initialise(options,R.getbeta(),A.getLoci(),PopulationLabels);

  options->PrintOptions();//NB: call after Populations is set		
  /*-------------------------------------------------------
  |  single individual, one population, allele frequencies |
  ---------------------------------------------------------*/
  if( options->getAnalysisTypeIndicator() == -1 && options->getPopulations() == 1 && strlen(options->getAlleleFreqFilename()) )
    IC->getOnePopOneIndLogLikelihood(&Log,A.getLoci(),PopulationLabels);

  /*----------
  | OTHERWISE |
  -----------*/
  //Initialise test objects
  else{
    Scoretest.Initialise(options, IC, A.getLoci(),chrm,&Log);  
    StratTest.Initialize( options, *(A.getLoci()) ,&Log);
    DispTest.Initialise(options,&Log, A.GetNumberOfCompositeLoci());

    //Initialise Output Files	
    if( options->getTextIndicator() ){
      L.InitializeOutputFile(PopulationLabels);
      R.InitializeOutputFile(options, IC,PopulationLabels);
      A.InitializeOutputFile(options, PopulationLabels);
      InitializeErgodicAvgFile(options,IC, &Log,&avgstream,PopulationLabels);
      Scoretest.InitialiseAssocScoreFile(PopulationLabels);  
    }
	
    SumLogTheta.SetNumberOfElements( options->getPopulations());
    IC->PreUpdate(L.getrhoalpha(),L.getrhobeta(),options);
 /*------------
  |  MAIN LOOP |
  ------------*/
    for( int iteration = 0; iteration <= options->getTotalSamples(); iteration++ ){
//Resets before updates
      Scoretest.Reset();
      A.Reset();
      SumLogTheta.SetElements( 0.0 );
	
//Updates  
      IC->Update(iteration,&SumLogTheta, R.getlambda(), R.getNoCovariates(), R.getbeta(),poptheta, options,
		 LociCorrSummary, A.getLoci(), chrm, L.getalpha(), _symmetric, _admixed, L.getrhoalpha(), L.getrhobeta(),
		 &LogFileStream, &MargLikelihood);

       A.UpdateAlleleFreqs(iteration,options->getBurnIn());
     if( iteration > options->getBurnIn() ){
       DispTest.UpdateBayesianPValueTest(*(A.getLoci()));
       if( options->getStratificationTest() )StratTest.calculate(IC, *(A.getLoci()));
     }  
     // Latent update should not need to take LociCorrSummary as an argument
     // Latent should not need to know anything about the number or positions of loci
     // with a global rho model, update of rho should be via a Metropolis random walk conditioned on the HMM likelihood
     // with a hierarchical rho model, update of hyperparameters should be via sufficient statistics: 
     // sum of rho and rho-squared over all individuals or gametes 
     L.Update(iteration, chrm, IC,&LociCorrSummary,&SumLogTheta,&poptheta,&LogFileStream);
     // Regression Update method should not need to take AnalysisTypeIndicator as an argument once the R object is initialized
     R.Update(options->getAnalysisTypeIndicator(),IC);

     if( iteration == options->getBurnIn() && options->getTestForAllelicAssociation() ){
       Scoretest.Output2(L.getalpha0(), &LogFileStream);
     }

     // output every 'getSampleEvery()' iterations
     if( !(iteration % options->getSampleEvery()) ){
       if( options->getAnalysisTypeIndicator() >= 0/*!= -3*/ ){//No Pop. Par. output for single individuals
	 Log.StartOutput(iteration,options->getTotalSamples());
	 L.OutputParams(iteration, &LogFileStream);
	 R.Output(iteration,&LogFileStream,options,IC);
	 A.OutputEta(iteration, options, &LogFileStream);
	 Log.EndOutput(iteration);
       }
       else 
	 cout << iteration << endl; 
      }     
     //Output and scoretest updates after BurnIn     
     if( iteration > options->getBurnIn() ){
       Scoretest.Update(IC->getTargetType(0),IC->getExpectedY0(),R.getlambda0());
       R.SumParameters(options->getAnalysisTypeIndicator());

       if( !(iteration % options->getSampleEvery()) ){
	 // output individual and locus parameters every 'getSampleEvery()' iterations
	 if ( strlen( options->getIndAdmixtureFilename() ) ) IC->accept();
	 if(options->getOutputAlleleFreq())A.accept();
       }//end of 'every' output

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
	   R.OutputErgodicAvg(samples,options, IC,&avgstream);
	   A.OutputErgodicAvg(samples,options,&avgstream);
	   if( options->getAnalysisTypeIndicator() == -1 ){
	     IC->OutputErgodicAvg(samples,&MargLikelihood,&avgstream);
	   }
	   avgstream << endl;
	 }
	 //Test output
	 if( options->getTestForDispersion() )  DispTest.Output(iteration - options->getBurnIn(), *(A.getLoci()));
	 if( options->getStratificationTest() ) StratTest.Output();
	 Scoretest.Output(iteration,PopulationLabels);
       }//end of 'every'*10 output
     }//end if after BurnIn
   }//end main loop

   if( options->getAnalysisTypeIndicator() == -1 )MargLikelihood.Output(&LogFileStream);
   if( options->getAnalysisTypeIndicator() == -3 )IC->Output(&LogFileStream);
   Scoretest.ROutput();
  }//end else

  for(int i=0; i<chrm->size(); i++){
    delete (*chrm)(i);
  }
 
  delete IC;//must call explicitly as IndAdmixOutputter destructor finishes writing to indadmixture.txt

//causes a crash
//    delete[] PopulationLabels;

  ProcessingTime(&Log, StartTime);
}

int main( int argc , char** argv ){
    PrintCopyrightNotice();

    int    xargc = argc;
    char **xargv = argv;    

    if (argc < 2) {
        cout << "Please specify an options file or command-line arguments\n";
        exit(1); // **need to print a warning message or list of user options to screen here
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
    str.erase(str.find_first_of(" \t\n\r"),str.find("=")-str.find_first_of(" \t\n\r"));//before '='
    str.erase(str.find_first_of(" \t\n\r"),str.find_last_of(" \t\n")-str.find_first_of(" \t\n\r")+1);//after '='
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
  long EndTime = time(0);
  tm timer;
  timer = *localtime( &EndTime );

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

  EndTime -= StartTime;
  timer = *localtime( &EndTime );
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
}

