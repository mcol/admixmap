#include "AlleleFreqs.h"
#include "StringSplitter.h"

//belongs in InputData
static void getLabels(const Vector_s& data, Vector_i temporary, string *labels)
{
    for (size_t i = 0, index = 0; i < data.size(); ++i) {
        if (temporary.GetNumberOfElements() == 1 || temporary(i)) {            
            labels[index++] = data[i];
        }
    }
}

AlleleFreqs::AlleleFreqs(){
  Vector_d null_Vector_d(1);
  Vector_i null_Vector_i(1);

  eta = null_Vector_d;
  etastep = null_Vector_d;
  etastep0 = 0.0;
  psi = null_Vector_d;
  tau = null_Vector_d; 
  SumEta = null_Vector_d;
  SumAcceptanceProb = null_Vector_d; 
  psi0 = 0.0;
  TuneEtaSampler = 0;
  w = 0;
  NumberAccepted = null_Vector_i;
  Number  = 0;
  Populations = 0;
  allelefreqoutput = 0;
  
}

AlleleFreqs::~AlleleFreqs(){

  for(int i=0; i<Loci.GetNumberOfCompositeLoci(); i++){
    delete Loci(i);
  }

  delete allelefreqoutput;

  if( isHistoricAlleleFreq ){
    delete [] TuneEtaSampler;
  }
}

void AlleleFreqs::Initialise(AdmixOptions *options,const Matrix_d& etaprior,LogWriter *Log,
			     std::string *PopulationLabels, double rho){
  Number = 0;
  Populations = options->getPopulations();
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ) isHistoricAlleleFreq = true;
  else isHistoricAlleleFreq = false;

  if( options->getOutputAlleleFreq() ){
    allelefreqoutput = new AlleleFreqOutputter(options,PopulationLabels);
  }

  psi.SetNumberOfElements( Populations );
  tau.SetNumberOfElements( Populations );

  SumAcceptanceProb.SetNumberOfElements( Populations );
  SumEta.SetNumberOfElements( Populations );

  LociCorrSummary.SetNumberOfElements( Loci.GetNumberOfCompositeLoci() );
  for( int j = 1; j < Loci.GetNumberOfCompositeLoci(); j++ )
    LociCorrSummary(j) = strangExp( -Loci.GetDistance( j ) * rho );

  // Matrix etaprior(1,1);
  Vector_d maxeta( Populations );
  if( isHistoricAlleleFreq ){
    w = 10;
    etastep0 = 2.0;
    etastep.SetNumberOfElements( Populations );
    etastep.SetElements( etastep0 );
    eta.SetNumberOfElements( Populations );
    NumberAccepted.SetNumberOfElements( Populations );
    TuneEtaSampler = new TuneRW[ Populations ];
    for( int k = 0; k < Populations; k++ )
      TuneEtaSampler[k].SetParameters( w, etastep0, 0.1, 100, 0.44 );
    if( strlen(options->getEtaPriorFilename()) ){
      Log->logmsg(true,"Loading gamma prior parameters for allele frequency dispersion from ");
      Log->logmsg(true,options->getEtaPriorFilename());
      Log->logmsg(true,".\n");
      //const Matrix_d& etaprior = data_->getEtaPriorMatrix();

      for( int k = 0; k < Populations; k++ ){
	psi(k) = etaprior( k, 0 );
	tau(k) = etaprior( k, 1 );
	Log->logmsg(true, "Population ");
	Log->logmsg(true, k);
	Log->logmsg(true, ": ");
	Log->logmsg(true, psi(k));
	Log->logmsg(true, " ");
	Log->logmsg(true, tau(k));
	Log->logmsg(true, "\n");
      }
    }
    else{
      psi.SetElements( 2 ); // default gamma prior with mean 400 
      tau.SetElements( 0.005 );
    }
    for( int k = 0; k < Populations; k++ ){
      for( int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
	maxeta(k) =  Loci(j)->GetPriorAlleleFreqs( k ).Sum();
	if( maxeta(k) > eta(k) ){
	  eta(k) = maxeta(k);
	}
      }
    }
  }

  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
  //Open output file for eta
  if ( strlen( options->getEtaOutputFilename() ) ){
    outputstream.open( options->getEtaOutputFilename(), ios::out );
    if( !outputstream )
      {
	Log->logmsg(true,"ERROR: Couldn't open dispparamfile\n");
	//exit( 1 );
      }
    else{
      Log->logmsg(true,"Writing dispersion parameter to ");
      Log->logmsg(true,options->getEtaOutputFilename());
      Log->logmsg(true,"\n");
    }
  }
  else{
    Log->logmsg(true,"No dispparamfile given\n");
    //exit(1);
  }

  }
  OpenFSTFile(options,Log);
  Log->logmsg(true,"Effective length of autosomes under study: ");
  Log->logmsg(true,Loci.GetLengthOfGenome());
  Log->logmsg(true," Morgans.\n");

  if( Loci.isX_data() ){
    Log->logmsg(true,"Effective length of X chromosome under study: ");
    Log->logmsg(true, Loci.GetLengthOfXchrm());
    Log->logmsg(true," Morgans.\n");
   }
}

void AlleleFreqs::load_f(double rho,Genome *chrm){
  int locus = 0;
  for( int j = 0; j < chrm->size(); j++ ){
    locus++;
    for( int jj = 1; jj < (*chrm)(j)->GetSize(); jj++ ){
      LociCorrSummary(locus) = exp( -Loci.GetDistance( locus ) * rho );
      locus++;
    }
  }
}

// Method samples allele frequency and prior allele frequency
// parameters.
void AlleleFreqs::UpdateAlleleFreqs(int iteration,int BurnIn){
  if( Loci(0)->IsRandom() ){
     
    Vector_d EtaParameters(3), probs;
    Matrix_d stats;
    MatrixArray_d data( Loci.GetNumberOfCompositeLoci() );
    
    double etanew, LogPostRatio;
    
    // Sample for prior frequency parameters
    if(isHistoricAlleleFreq ){
      Loci.SamplePriorAlleleFreqs( eta );
    }
    // Sample for allele frequencies
    Loci.SampleAlleleFreqs( 1 );
    
    // Sample for allele frequency dispersion parameters, eta, using
    // Metropolis random-walk.
    if(  isHistoricAlleleFreq ){
      Number++;
      for( int k = 0; k < Populations; k++ ){
	// Sample eta from truncated log-normal distribution.
	vector< Vector_d > munew;
	do{
	  etanew = exp( gennor( log( eta(k) ), etastep(k) ) );
	}while( etanew > 5000.0 );
	// Prior log-odds ratio         
	LogPostRatio = ( psi(k) - 1 ) * (log(etanew) - log(eta(k)))
	  - tau(k) * ( etanew - eta(k) );
	// Log-likelihood ratio; numerator of integrating constant
	LogPostRatio += 2 * Loci.GetNumberOfCompositeLoci()
	  * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( eta(k) ) );
	for( int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
	  Vector_d mu = Loci(j)->GetPriorAlleleFreqs( k );
	  munew.push_back( mu * etanew / eta(k) );
	  Vector_d SumLogFreqs = Loci(j)->GetStatsForEta( k );
	  for( int l = 0; l < Loci(j)->GetNumberOfStates(); l++ ){
	    // Denominator of integrating constant
	    LogPostRatio += 2*(gsl_sf_lngamma( mu(l) ) - gsl_sf_lngamma( munew[j](l) ));
	    // SumLogFreqs = log phi_1 + log phi_2
	    LogPostRatio += (munew[j](l) - mu(l))*SumLogFreqs(l);
	  }
	}
	
	// Log acceptance probability = Log posterior ratio since the
	// proposal ratio (log-normal) cancels with prior.
	
	// Acceptance test.
	if( log( myrand() ) < LogPostRatio ){
	  eta(k) = etanew;
	  Loci.UpdatePriorAlleleFreqsGlobal( k, munew );
	  SumAcceptanceProb(k)++;
	  NumberAccepted(k)++;
	}
	
	if( !( Number % w ) ){
	  etastep(k) = TuneEtaSampler[k].UpdateSigma( NumberAccepted(k) );
	  NumberAccepted(k) = 0;
	}
      }
      
      if( !( Number % w ) ){
	Number = 0;
      }
      
      if( iteration > BurnIn )
	SumEta += eta;
    }
    
    if( iteration > BurnIn && isHistoricAlleleFreq ){
      Loci.UpdateFst();
    }
  }
}

void AlleleFreqs::InitializeOutputFile(AdmixOptions *options, std::string *PopulationLabels)
{
  // Header line of paramfile
  if( options->getAnalysisTypeIndicator() >= 0 ){

    //Dispersion parameters (eta)
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      for( int k = 0; k < Populations; k++ ){
	outputstream << "\"eta." << PopulationLabels[k].substr(1);
      }
    }
    outputstream << endl;
  }
}

void AlleleFreqs::OutputErgodicAvg( int samples,AdmixOptions *options, std::ofstream *avgstream)
{
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    for( int j = 0; j < Populations; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumEta(j) / samples << " ";
    }
  }
}

void AlleleFreqs::OutputEta(int iteration, AdmixOptions *options, std::ofstream *LogFileStreamPtr){
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
  //output to logfile
    if( !options->useCOUT() || iteration == 0 )
      {
	for( int j = 0; j < Populations; j++ ){
	  LogFileStreamPtr->width(9);
	  (*LogFileStreamPtr) << setprecision(6) << eta(j) << " ";
	}
      }
  //output to screen
    if( options->useCOUT() )
      {
	for( int j = 0; j < Populations; j++ ){
	  cout.width(9);
	  cout << setprecision(6) << eta(j) << " ";
	}
      }
  //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      for( int j = 0; j < Populations; j++ ){
	outputstream.width(9);
	outputstream << setprecision(6) << eta(j);
      }
     outputstream << endl;
    }
  }
}

Vector_d *AlleleFreqs::geteta(){
  return &eta;
}
Vector_d *AlleleFreqs::getSumEta(){
  return &SumEta;
}
Genome *AlleleFreqs::getLoci(){
  return &Loci;
}
int AlleleFreqs::GetNumberOfCompositeLoci(){
  return Loci.GetNumberOfCompositeLoci();
}
void AlleleFreqs::accept(){
if( Loci(0)->IsRandom() ){
	   Loci.accept(*allelefreqoutput);
	 }
}
void AlleleFreqs::OutputFST(bool IsPedFile){
  for( int j = 0; j < GetNumberOfCompositeLoci(); j++ ){
    if(IsPedFile)
      fstoutputstream << "\"" << Loci(j)->GetLabel(0) << "\"";
    else
      fstoutputstream << Loci(j)->GetLabel(0);
    fstoutputstream << " " << Loci(j)->GetFst() << endl;
  }
}

void AlleleFreqs::Reset(){
  Loci.ResetScoreForMisSpecOfAlleleFreqs();
  Loci.ResetLikelihoodAlleleFreqs();
}
void AlleleFreqs::OpenFSTFile(AdmixOptions *options,LogWriter *Log){
  if( options->getOutputFST() ){
    Log->logmsg(true,options->getFSTOutputFilename());
    Log->logmsg(true,"\n");
    fstoutputstream.open( options->getFSTOutputFilename(), ios::out );
    if( !fstoutputstream ){
      Log->logmsg(true,"ERROR: Couldn't open fstoutputfile\n");
      exit( 1 );
    }
  }
}
void AlleleFreqs::loadAlleleStatesAndDistances(vector<string> * ChrmLabels,AdmixOptions *options,InputData *data_, LogWriter *Log){
  string *LociLabelsCheck = 0;

  // Load number of allelic states and distances.
  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getGeneInfoFilename());
  Log->logmsg(false,".\n");

  Matrix_d& locifileData = (Matrix_d&) data_->getGeneInfoMatrix();
  
  LociLabelsCheck = new string[ locifileData.GetNumberOfRows() ];
  int numCompLoci = 0;
  for( int i = 0; i < locifileData.GetNumberOfRows(); i++ )
    if( locifileData( i, 1 ) == 0.0 )
      numCompLoci++;
  Loci.SetNumberOfCompositeLoci(locifileData.GetNumberOfRows() - numCompLoci);
  for(int i=0;i < GetNumberOfCompositeLoci();i++){
    Loci(i) = new CompositeLocus();
  }

  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getGeneticDataFilename());
  Log->logmsg(false,".\n");


  Loci.SetNumberOfLoci(1);

  // Set number of alleles at each locus
  int index =0;
  size_t next_line = 0;
  for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
    ++next_line;

    const Vector_s& m = data_->getGeneInfoData()[next_line];

    LociLabelsCheck[index] = m[0];
    if (m.size() == 4)
      ChrmLabels->push_back(m[3]);
    Loci(i)->SetNumberOfAllelesOfLocus( 0, (int)locifileData( i, 0 ) );
    Loci.SetDistance( i, locifileData( index, 1 ) );
    while( index < locifileData.GetNumberOfRows() - 1 && locifileData( index + 1, 1 ) == 0 ){
      ++next_line;

      Loci(i)->AddLocus( (int)locifileData( index + 1, 0 ) );
      index++;
      LociLabelsCheck[index] = m[0];
    }
    CompositeLocus *locus = (CompositeLocus*)Loci(i);
    locus->SetNumberOfLabels();
    index++;
    //Log->logmsg(false,Loci(i)->GetNumberOfLoci());
    //Log->logmsg(false," ");
  }
  Log->logmsg(false,"\n");

  if( options->getTextIndicator() ){

    Vector_s labels = data_->getGeneticData()[0];

    Vector_d vtemp = locifileData.GetColumn(1);
    Log->logmsg(true, vtemp.GetNumberOfElements());Log->logmsg(true," simple loci\n");
    vtemp.AddElement(0); // Forces SetLabels method to ignore first row of loci.txt (GenotypesFile)
    // Add a sex column if it is not included
    if( ! options->genotypesSexColumn() ){
      labels.insert(labels.begin(), "\"extracol\"");
    }
    Loci.SetLabels(labels, vtemp);
  }

  index = 0;
  for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
    for( int j = 0; j < Loci(i)->GetNumberOfLoci(); j++ ){
      if( LociLabelsCheck[index].compare( Loci(i)->GetLabel(0) ) ){
	Log->logmsg(true, "Error in loci names in genotypes file and loci file at loci\n" );
	Log->logmsg(true, i);
	Log->logmsg(true, LociLabelsCheck[index] );
	Log->logmsg(true, " " );
	Log->logmsg(true, j);
	Log->logmsg(true, Loci(i)->GetLabel(j) );
	Log->logmsg(true, "\n" );
	exit(0);
      }
      index++;
    }
  }  
  delete [] LociLabelsCheck;
}

void AlleleFreqs::LoadAlleleFreqs(AdmixOptions *options, Genome **chrm,LogWriter *Log, InputData *data_,std::string **PopulationLabels)
{
  int newrow;
  int row = 0;

  Matrix temp2, elements;
  Matrix_d temporary;
  vector<string> ChrmLabels;

  Populations = options->getPopulations();
  checkLociNames(options,data_);
  loadAlleleStatesAndDistances(&ChrmLabels,options,data_, Log);

  // Load allele frequencies
  if( strlen( options->getAlleleFreqFilename() ) ){

    Log->logmsg(false,"Loading ");
    Log->logmsg(false,options->getAlleleFreqFilename());
    Log->logmsg(false,".\n");

    temporary = data_->getAlleleFreqMatrix();

    if(temporary.GetNumberOfRows()-1 != Loci.GetNumberOfStates()-GetNumberOfCompositeLoci()){
      Log->logmsg(true,"Incorrect number of rows in allelefreqsfile.\n");
      Log->logmsg(true,"Expecting ");
      Log->logmsg(true,Loci.GetNumberOfStates()-GetNumberOfCompositeLoci()+1);
      Log->logmsg(true," rows, where as there are ");
      Log->logmsg(true,temporary.GetNumberOfRows());
      Log->logmsg(true," rows.\n");
      exit(0);
    }
    //options->setPopulations( temporary.GetNumberOfCols() - options->getTextIndicator() );
    Populations = temporary.GetNumberOfCols() - options->getTextIndicator();
    if( options->getTextIndicator() ){
      temporary = temporary.SubMatrix( 1, temporary.GetNumberOfRows() - 1, 1, Populations );
 
      *PopulationLabels = new string[ Populations ];

      Vector_i vtemp( Populations + 1 );
      vtemp.SetElements( 1 );
      vtemp(0) = 0;
      ::getLabels(data_->getAlleleFreqData()[0], vtemp, *PopulationLabels);
    }

    for( int i = 0; i < GetNumberOfCompositeLoci(); i++ )
      {
	newrow = row + Loci(i)->GetNumberOfStates() - 1;
	Loci(i)->SetAlleleFreqs( (temporary.Double()).SubMatrix( row, newrow - 1, 0, Populations - 1 ) );
	row = newrow;
      }

    if( row != temporary.GetNumberOfRows() ){
      Log->logmsg(true,"Inconsistency in ");
      Log->logmsg(true,options->getAlleleFreqFilename());
      Log->logmsg(true," and ");
      Log->logmsg(true,options->getGeneInfoFilename());
      Log->logmsg(true,"\n");
      Log->logmsg(true,row);
      Log->logmsg(true," ");
      Log->logmsg(true,temporary.GetNumberOfRows());
      Log->logmsg(true,"\n");
      exit(0);
    }
  }
  else if( strlen( options->getHistoricalAlleleFreqFilename() ) || strlen( options->getPriorAlleleFreqFilename() ) ){
    const Vector_s* alleleFreqLabels = 0;
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      alleleFreqLabels = &data_->getHistoricalAlleleFreqData()[0];
      Log->logmsg(false,"Loading ");
      Log->logmsg(false,options->getHistoricalAlleleFreqFilename());
      Log->logmsg(false,".\n");
      temporary = data_->getHistoricalAlleleFreqMatrix();
      //options->setPopulations(temporary.GetNumberOfCols() - options->getTextIndicator());
       Populations = temporary.GetNumberOfCols() - options->getTextIndicator();
       
      if( temporary.GetNumberOfRows() != Loci.GetNumberOfStates()+1 ){
	Log->logmsg(true,"Incorrect number of rows in historicalallelefreqsfile.\n");
	Log->logmsg(true,"Expecting ");
	Log->logmsg(true,Loci.GetNumberOfStates()+1);
	Log->logmsg(true," rows, but there are ");
	Log->logmsg(true,temporary.GetNumberOfRows());
	Log->logmsg(true," rows.\n");
	exit(0);
      }
    } else {
      alleleFreqLabels = &data_->getPriorAlleleFreqData()[0];
      Log->logmsg(false,"Loading ");
      Log->logmsg(false,options->getPriorAlleleFreqFilename());
      Log->logmsg(false,".\n");
      temporary = data_->getPriorAlleleFreqMatrix();
      //options->setPopulations(temporary.GetNumberOfCols() - options->getTextIndicator());
       Populations = temporary.GetNumberOfCols() - options->getTextIndicator();

      if( temporary.GetNumberOfRows() != Loci.GetNumberOfStates()+1 ){
	Log->logmsg(true,"Incorrect number of rows in priorallelefreqsfile.\n");
	Log->logmsg(true,"Expecting ");
	Log->logmsg(true,Loci.GetNumberOfStates()+1);
	Log->logmsg(true," rows, where as there are ");
	Log->logmsg(true,temporary.GetNumberOfRows());
	Log->logmsg(true," rows.\n");
	exit(1);
      }
    }
 
    temporary = temporary.SubMatrix( 1, temporary.GetNumberOfRows() - 1, 1, Populations );
    *PopulationLabels = new string[ Populations ];

    Vector_i vtemp( Populations + 1 );
    vtemp.SetElements( 1 );
    vtemp(0) = 0;
    ::getLabels(*alleleFreqLabels, vtemp, *PopulationLabels);

    for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
      newrow = row + Loci(i)->GetNumberOfStates();
      if( strlen( options->getHistoricalAlleleFreqFilename() ) )
	Loci(i)->SetHistoricalAlleleFreqs( temporary.SubMatrix( row, newrow - 1, 0,
								Populations - 1 ) );
      else
	Loci(i)->SetPriorAlleleFreqs( temporary.SubMatrix( row, newrow - 1, 0, Populations - 1 ), options->getFixedAlleleFreqs() );
      row = newrow;
    }
  }
  else{
    Loci.SetDefaultAlleleFreqs( Populations );
    *PopulationLabels = new string[ Populations ];
    for( int j = 0; j < Populations; j++ ){
      stringstream poplabel;
      string result;
      poplabel << "\"A" << j << "\"";
      result = poplabel.str();
      (*PopulationLabels)[j] = result;
    }
      
  }
 *chrm = Loci.GetChromosomes( Populations, ChrmLabels );
 options->setPopulations(Populations);
 //(**)
  Log->logmsg(false,Loci.size());
  Log->logmsg(false," loci; ");
  Log->logmsg(false,(*chrm)->size());
  Log->logmsg(false," chromosomes\n");
}

void AlleleFreqs::checkLociNames(AdmixOptions *options,InputData *data_){
  // Check that loci labels in locusfile are unique and that they match the names in the genotypes file.
  
    const Matrix_s& geneInfoData = data_->getGeneInfoData();;
    const Matrix_s& geneticData  = data_->getGeneticData();

    // Check loci names uniqueness.    
    for (size_t i = 1; i < geneInfoData.size(); ++i) {
        for (size_t j = i + 1; j < geneInfoData.size(); ++j) {   
            if (geneInfoData[i][0] == geneInfoData[j][0]) {
                    cerr << "Error in locusfile. Two different loci have the same name. "
                         << geneInfoData[i][0] << endl;
                    exit(2);            
            }
        }
    }

    const size_t numLoci = geneInfoData.size() - 1;

    // Determine if "Sex" column present in genotypes file.
    if (numLoci == geneticData[0].size() - 1) {
        options->genotypesSexColumn(0);
    } else if (numLoci == geneticData[0].size() - 2) {
        options->genotypesSexColumn(1);
    } else {
        cerr << "Error. Number of loci in genotypes file does not match number in locus file." << endl;
        exit(2);
    }

    // Compare loci names in locus file and genotypes file.
    for (size_t i = 1; i <= numLoci; ++i) {
        if (geneInfoData[i][0] != geneticData[0][i + options->genotypesSexColumn()]) {
            cout << "Error. Loci names in locus file and genotypes file are not the same." << endl;
            cout << "Loci names causing an error are: " << geneInfoData[i][0] << " and " 
                 << geneticData[0][i + options->genotypesSexColumn()] << endl;
            cout << options->genotypesSexColumn() << endl;
            exit(2);
        }
    } 
}

void AlleleFreqs::getLabels( const string buffer, Vector_i temporary, string *labels )
{
    StringSplitter splitter;
    const Vector_s& labels_tmp = splitter.split(buffer);

    for (size_t i = 0, index = 0; i < labels_tmp.size(); ++i) {
        if (temporary.GetNumberOfElements() == 1 || temporary(i)) {            
            labels[index++] = labels_tmp[i];
        }
    }
}

Vector_d AlleleFreqs::getLociCorrSummary(){
  return LociCorrSummary;
}
//generic method - should be elsewhere?
double AlleleFreqs::strangExp( double x )
{
  double y;
  if( x > -700 )
    y = exp(x);
  else
    y = 0;
  return( y );
}
