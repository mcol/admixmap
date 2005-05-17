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
  eta = 0;
  etastep = 0;
  etastep0 = 0.0;
  psi = 0;
  tau = 0; 
  SumEta = 0;
  SumAcceptanceProb = 0; 
  psi0 = 0.0;
  TuneEtaSampler = 0;
  w = 0;
  NumberAccepted = 0;
  Number  = 0;
  Populations = 0;
  allelefreqoutput = 0;
  RandomAlleleFreqs = 0;
  IsHistoricAlleleFreq = false;
  AlleleCounts = 0;
  Freqs = 0;
  AlleleFreqsMAP = 0;
  HistoricAlleleFreqs = 0;
  HistoricLikelihoodAlleleFreqs = 0;
  PriorAlleleFreqs = 0;
  SumAlleleFreqs = 0;
  Fst = 0;
  SumFst = 0;
  allelefreqoutput = 0;
  TuneEtaSampler = 0;
  LociCorrSummary = 0;
  pp = 0;  
}

AlleleFreqs::~AlleleFreqs(){

  for(unsigned int i=0; i < Loci.GetNumberOfCompositeLoci(); i++){
    delete Loci(i);
  }
  delete allelefreqoutput;

  delete[] TuneEtaSampler;

  delete[] Freqs;
  delete[] AlleleCounts;
  delete[] PriorAlleleFreqs;
  delete[] HistoricLikelihoodAlleleFreqs;
  delete[] SumAlleleFreqs;
  delete[] AlleleFreqsMAP;
  delete[] HistoricAlleleFreqs;
  delete[] Fst;
  delete[] SumFst;
  delete[] psi;
  delete[] tau;
  delete[] NumberAccepted;
  delete[] SumAcceptanceProb;
  delete[] etastep;
  delete[] SumEta;
  delete[] LociCorrSummary;
  delete[] pp;
}

void AlleleFreqs::loadAlleleStatesAndDistances(vector<string> * ChrmLabels,AdmixOptions *options,InputData *data_, LogWriter *Log){
  string *LociLabelsCheck = 0;

  // Load number of allelic states and distances.
  //also sets number of composite loci and initialises the allelefreq arrays
  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getGeneInfoFilename());
  Log->logmsg(false,".\n");

  Matrix_d& locifileData = (Matrix_d&) data_->getGeneInfoMatrix();
  
  LociLabelsCheck = new string[ locifileData.GetNumberOfRows() ];
  int numCompLoci = 0;
  // numCompLoci counts the number of simple Loci that are in neither simple loci nor the first in a compound locus
  // would be neater to do some thing like
  //numCompLoci = locifileData.GetNumberOfRows()
  //  for( int i = 0; i < locifileData.GetNumberOfRows(); i++ )
  //   if( locifileData( i, 1 ) == 0.0 ) numCompLoci--;
  //so that numCompLoci is the number of Composite Loci
  // could also store this at class scope to save calls to Loci::GetNumberOfCompositeLoci

  for( int i = 0; i < locifileData.GetNumberOfRows(); i++ )
    if( locifileData( i, 1 ) == 0.0 ) numCompLoci++;

  //determine number of composite loci and allocate arrays
  // should only allocate the arrays that are needed
  //eg we don't always use HistoricAlleleFreqs
  int NCL = locifileData.GetNumberOfRows() - numCompLoci;
  Loci.SetNumberOfCompositeLoci(NCL);
  Freqs = new Matrix_d[NCL];
  AlleleFreqsMAP = new Matrix_d[NCL];
  HistoricAlleleFreqs = new Matrix_d[NCL];
  AlleleCounts = new Matrix_i[NCL];
  HistoricLikelihoodAlleleFreqs = new Matrix_d[NCL];
  PriorAlleleFreqs = new Matrix_d[NCL];
  SumAlleleFreqs = new Matrix_d[NCL];
  MuProposal = new std::vector<TuneRW>[NCL];


  for(int i=0;i < GetNumberOfCompositeLoci();i++){
    Loci(i) = new CompositeLocus();
    Loci(i)->SetNumberOfLoci(1);
  }

  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getGeneticDataFilename());
  Log->logmsg(false,".\n");

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

  // checks of input data files should be in class InputData
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

void AlleleFreqs::LoadAlleleFreqs(AdmixOptions *options, Chromosome ***chrm,LogWriter *Log, InputData *data_,std::string **PopulationLabels)
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

    InitialiseAlleleFreqs( (temporary.Double()), Populations);
//     for( int i = 0; i < GetNumberOfCompositeLoci(); i++ )
//       {
// 	newrow = row + Loci(i)->GetNumberOfStates() - 1;
// 	Loci(i)->SetAlleleFreqs( (temporary.Double()).SubMatrix( row, newrow - 1, 0, Populations - 1 ) );
// 	row = newrow;
//       }

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
	InitialiseHistoricAlleleFreqs( temporary.SubMatrix( row, newrow - 1, 0,Populations - 1 ), i );
      else
	InitialisePriorAlleleFreqs( temporary.SubMatrix( row, newrow - 1, 0, Populations - 1 ), i,options->getFixedAlleleFreqs());
      row = newrow;
    }
  }
  else{
    //Loci.SetDefaultAlleleFreqs( Populations );
    SetDefaultAlleleFreqs( Populations );
    *PopulationLabels = new string[ Populations ];
    for( int j = 0; j < Populations; j++ ){
      stringstream poplabel;
      string result;
      poplabel << "\"A" << j << "\"";
      result = poplabel.str();
      (*PopulationLabels)[j] = result;
    }
      
  }

  (*chrm) = Loci.GetChromosomes(Populations, ChrmLabels );
  
  Loci.SetSizes();
  options->setPopulations(Populations);
  pp = new double[Populations];
  //(**)
  Log->logmsg(false,Loci.GetNumberOfCompositeLoci());
  Log->logmsg(false," loci; ");
  Log->logmsg(false, Loci.GetNumberOfChromosomes());
  Log->logmsg(false," chromosomes\n");
}

void AlleleFreqs::checkLociNames(AdmixOptions *options,InputData *data_){
  // Check that loci labels in locusfile are unique and that they match the names in the genotypes file.
  
    const Matrix_s& geneInfoData = data_->getGeneInfoData();;
    const Matrix_s& geneticData  = data_->getGeneticData();

    // Check loci names are unique    
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

void AlleleFreqs::Initialise(AdmixOptions *options,const Matrix_d& etaprior,LogWriter *Log,
			     std::string *PopulationLabels, double rho){
  Number = 0;
  Populations = options->getPopulations();
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ) IsHistoricAlleleFreq = true;
  else IsHistoricAlleleFreq = false;

  if(IsRandom() &&  options->getOutputAlleleFreq() ){
    allelefreqoutput = new AlleleFreqOutputter(options,PopulationLabels);
  }

  //set up alleleprobs and hap pair probs
  //NB: HaplotypePairProbs in Individual must be set first
  for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
    Loci(i)->Initialise(Freqs[i]);
  }

  // these lines should be moved to Genome class ?
  LociCorrSummary = new double[ Loci.GetNumberOfCompositeLoci() ];
  for(unsigned int j = 1; j < Loci.GetNumberOfCompositeLoci(); j++ )
     LociCorrSummary[j] = ( -Loci.GetDistance( j ) * rho > -700) ? exp( -Loci.GetDistance( j ) * rho ) : 0.0;

  // settings for sampling of dispersion parameter
  // Matrix etaprior(1,1);
  psi = new double[ Populations ];
  tau = new double[ Populations ];

  SumAcceptanceProb = new double[ Populations ];
  SumEta = new double[ Populations ];

  Vector_d maxeta( Populations );
  if( IsHistoricAlleleFreq ){
    w = 10;
    etastep0 = 2.0;
    etastep =new double[ Populations ];
    for(int k=0;k<Populations;++k)etastep[k] = etastep0;

    eta = new double[Populations];
    NumberAccepted =new int[ Populations ];
    TuneEtaSampler = new TuneRW[ Populations ];
    for( int k = 0; k < Populations; k++ )
      TuneEtaSampler[k].SetParameters( w, etastep0, 0.1, 100, 0.44 );
    if( strlen(options->getEtaPriorFilename()) ){
      Log->logmsg(true,"Loading gamma prior parameters for allele frequency dispersion from ");
      Log->logmsg(true,options->getEtaPriorFilename());
      Log->logmsg(true,".\n");
      //const Matrix_d& etaprior = data_->getEtaPriorMatrix();

      for( int k = 0; k < Populations; k++ ){
	psi[k] = etaprior( k, 0 );
	tau[k] = etaprior( k, 1 );
	Log->logmsg(true, "Population ");
	Log->logmsg(true, k);
	Log->logmsg(true, ": ");
	Log->logmsg(true, psi[k]);
	Log->logmsg(true, " ");
	Log->logmsg(true, tau[k]);
	Log->logmsg(true, "\n");
      }
    }
    else{
      //psi.SetElements( 2 ); // default gamma prior with mean 400 
      //tau.SetElements( 0.005 );
      for(int k=0; k<Populations; ++k){
	psi[k] = 2.0;
	tau[k] = 0.005;
      }
    }
    for( int k = 0; k < Populations; k++ ){
      
//       // old method; sets eta to the sum of priorallelefreqs
//       for( int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
//        	maxeta(k) =  GetPriorAlleleFreqs(j,k).Sum();
//        	if( maxeta(k) > eta[k] ){
//        	  eta[k] = maxeta(k);
//        	}
//       }
      
      //Initialise eta at its prior expectation
      eta[k] = psi[k]/tau[k];
      //Rescale priorallelefreqs so the columns sum to eta 
      for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ )
	PriorAlleleFreqs[j].SetColumn(k, PriorAlleleFreqs[j].GetColumn(k) * eta[k] / PriorAlleleFreqs[j].GetColumn(k).Sum());
     }
  
    //Open output file for eta
    if ( options->getIndAdmixHierIndicator() && strlen( options->getEtaOutputFilename() ) ){
      InitializeEtaOutputFile(options, PopulationLabels, Log); 
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

void AlleleFreqs::load_f(double rho,Chromosome **chrm){
  int locus = 0;
  for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
    locus++;
    for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
      LociCorrSummary[locus] = exp( -Loci.GetDistance( locus ) * rho );
      locus++;
    }
  }
}

void AlleleFreqs::InitialiseAlleleFreqs(Matrix_d NewAlleleFreqs, int Pops){
/**
 * Sets the frequencies of each haplotype at each composite locus.
 * obsolescent - maintained for compatibility with old format 
 * in which fixed allele freqs were specified in input data with last allele omitted 
 * 
 * NewAlleleFreqs - a matrix containing allele frequencies, read from input data
 * Rows index alleles or haplotypes, omitting the last allele
 * as its frequency is 1 minus sum of all other frequencies. 
 * Cols index populations.  Thus, for a composite locus with four states
 *   and European and African populations, the matrix might be:
 *
 *             Population
 *
 *            | EUR | AFR |
 *         ---|-----|-----|
 *          0 | 0.5 | 0.2 |
 *   State  1 | 0.1 | 0.2 |
 *          2 | 0.2 | 0.5 |
 * also initializes score test - this doesn't belong here
 */
  int newrow;
  int row = 0;
  int NumberOfStates;

  for( int i = 0; i < GetNumberOfCompositeLoci(); i++ )
    {
      NumberOfStates = Loci(i)->GetNumberOfStates();
     // check that number of alleles is correct 
      if( NewAlleleFreqs.GetNumberOfRows() != NumberOfStates - 1 ){//should use logmsg
	cout << "Error in number of alleles in SetAlleleFreqs.\n";
	cout << "Number of states = " << NumberOfStates << endl;
	cout << "AlleleFreqs has " << NewAlleleFreqs.GetNumberOfRows() << " rows.\n";
	exit(0);
      }
      else{
	Loci(i)->SetNumberOfPopulations(Pops);
	newrow = row + Loci(i)->GetNumberOfStates() - 1;
	// initialize Freqs
	Freqs[i] = NewAlleleFreqs.SubMatrix( row, newrow - 1, 0, Pops - 1 );
	// set size of allele counts matrix at this locus
	AlleleCounts[i].SetNumberOfElements(Loci(i)->GetNumberOfStates(),Pops);
	// initialize score test
	Loci(i)->InitialiseScoreTest(Pops); // ?move this to ScoreTests class
	row = newrow;

      }//end else
    }

}


void AlleleFreqs::InitialisePriorAlleleFreqs(Matrix_d New, int i, bool fixed){
/**
 *
 * Initialises the frequencies of each allele at in the ith
 * composite locus, given Dirichlet priors in matrix New.  Allele freqs
 * are set to their prior expectation 
 * If fixed, allele freqs will be fixed at their prior expectations   
 *
 * NewPriorAlleleFreqs - a matrix containing 
 *   parameters for the Dirichlet prior distribution of the allele frequencies. The first dimension is the allele number, 
 *   being in the range of zero to two less than the number of states
 *   [see GetNumberOfStates()]. The sum of the prior parameters over all alleles in a population 
 *   (sumalpha) can be interpreted as 
 +   the "prior sample size". The second dimension is the population. Thus, for a 
 *   composite locus with four states and European and African 
 *   populations, the matrix might be:
 *
 *             Population
 *
 *            | EUR | AFR |
 *         ---|-----|-----|
 *          0 | 9.0 | 3.0 |
 *   State  1 | 3.0 | 4.0 |
 *          2 | 1.0 | 8.0 |
 *          3 | 2.0 | 1.0 |
 */
  double sumalpha;
  int Pops = New.GetNumberOfCols();
  Loci(i)->SetNumberOfPopulations(Pops);
  // check size of prior allele freqs matrix
  if( New.GetNumberOfRows() != Loci(i)->GetNumberOfStates() ){
    cout << "Error in number of alleles in InitialisePriorAlleleFreqs.\n";
    cout << "Number of states " << Loci(i)->GetNumberOfStates() << endl;
    cout << "PriorAlleleFreqs has " << New.GetNumberOfRows() << " rows.\n";
  }
  else{
    // initialize Freqs and AlleleCounts
    // set size of allele freqs matrix for this locus
    Freqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates()-1, Pops);
    // set size of allele counts matrix
    AlleleCounts[i].SetNumberOfElements(Loci(i)->GetNumberOfStates(), Pops);

    // allele frequencies are initialised as expectations over the Dirichlet prior distribution, 
    // by dividing each prior parameter by the sum of the parameters.     
    for( int j = 0; j < Pops; j++ ){
      sumalpha = ( New.GetColumn(j) ).Sum();
      for( int k = 0; k < Loci(i)->GetNumberOfStates() - 1; k++ )
	Freqs[i]( k, j ) = ( New( k, j ) ) / sumalpha;
    }
    if(fixed){
      Loci(i)->InitialiseScoreTest(Pops);
    }
    else{
      PriorAlleleFreqs[i] = New;
      SumAlleleFreqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() -1, Pops);
      RandomAlleleFreqs = 1;
    }
  }
}

void AlleleFreqs::InitialiseHistoricAlleleFreqs(Matrix_d New, int i){
  /**
   * This method sets "historical allele frequencies", where the model has been specified to allow the 
   * allele freqs in the admixed population 
   * to vary from the historical allele frequencies in the unadmixed ancestral populations that have 
   * been sampled. 
   * Otherwise as for InitalisePriorAlleleFreqs
   * 
   */
  // unnecessary duplication of code - should use same method as above.  only difference is in 
  // how sumalpha is set
  double sumalpha;

  if( New.GetNumberOfRows() != Loci(i)->GetNumberOfStates() ){
    cout << "Error in number of alleles in SetHistoricalAlleleFreqs.\n";
    cout << "Number of states " << Loci(i)->GetNumberOfStates() << endl;
    cout << "HistoricalAlleleFreqs has "<< New.GetNumberOfRows() << " rows.\n";
  }
  HistoricLikelihoodAlleleFreqs[i] = New;
  PriorAlleleFreqs[i] = New + 0.501;
  int Pops = New.GetNumberOfCols();
  Loci(i)->SetNumberOfPopulations(Pops);

  // initialize Freqs
  Freqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() - 1, Pops);
  HistoricAlleleFreqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() - 1, Pops);
  AlleleCounts[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() , Pops);
  //delete[] SumEta;
  //SumEta = new double[Pops];
  for( int j = 0; j < Pops; j++ ){
    sumalpha = ( HistoricLikelihoodAlleleFreqs[i].GetColumn(j) ).Sum();
   //sumalpha = ( New.GetColumn(j) ).Sum();
    for( int k = 0; k < Loci(i)->GetNumberOfStates() - 1; k++ )
      Freqs[i]( k, j ) = ( New( k, j ) ) / sumalpha;
  }
  SumAlleleFreqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() - 1, Pops);
  RandomAlleleFreqs = 1;
  Fst = alloc2D_d(GetNumberOfCompositeLoci(), Pops);
  SumFst = alloc2D_d(GetNumberOfCompositeLoci(), Pops);
  // set size of vector MuProposal
  if( Loci(i)->GetNumberOfStates() > 2 ){
    MuProposal[i].resize( Populations );
    for( int k = 0; k < Populations; k++ ){
      MuProposal[i][k].SetParameters( 10, 0.01, 0.001, 0.1, 0.23 );
    }
  }
}

void AlleleFreqs::SetDefaultAlleleFreqs(int Pops){
/**
 * Given the number of ancestral populations, sets default values for
 * allele frequencies and prior allele frequencies.
 *
 * populations - the number of ancestral populations
 */
  // this check should be moved into AdmixOptions class
   if( Pops < 1 ){
      cout << "Error in SetDefaultAlleleFreqs( int populations ).\n";
      cout << "Number of populations = " << Pops << endl;
      exit(0);
   }
  // this check should be moved into InputData class
   for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
     if(Loci(i)->GetNumberOfStates() < 2){
       cout << "Error: The number of alleles at a locus is < 2. There must be at least two different alleles at each locus." << endl;
       exit(0);
     }
    // more duplicated code - should do this within method InitializePriorAlleleFreqs
     Loci(i)->SetNumberOfPopulations(Pops);
     PriorAlleleFreqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates(), Pops);
    // reference prior on allele freqs: all elements of parameter vector set to 0.5
    // this is unrealistic for large haplotypes - should set all elements to sum to 1
     PriorAlleleFreqs[i].SetElements(0.5);
     // initialize frequencies as equal for all alleles at locus
     Freqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() - 1, Pops);
     Freqs[i].SetElements(1.0/Loci(i)->GetNumberOfStates());

     AlleleCounts[i].SetNumberOfElements(Loci(i)->GetNumberOfStates(), Pops);
     SumAlleleFreqs[i].SetNumberOfElements(Loci(i)->GetNumberOfStates() - 1, Pops);
     RandomAlleleFreqs = 1;
   }
}

// Method samples allele frequency and prior allele frequency
// parameters.
void AlleleFreqs::Update(int iteration,int BurnIn){
  //Reset();

  if( IsRandom() ){
     
    Vector_d EtaParameters(3), probs;
    Matrix_d stats;
    MatrixArray_d data( Loci.GetNumberOfCompositeLoci() );
    
    double etanew, LogPostRatio;
    
    // Sample for prior frequency parameters mu, using eta, the sum of the frequency parameters for each locus.
    if(IsHistoricAlleleFreq ){
      for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
	if( Loci(i)->GetNumberOfStates() == 2 )
	  SamplePriorAlleleFreqs1D( i);
	else
	  SamplePriorAlleleFreqsMultiDim( i);
      }
    }

    // Sample allele frequencies and set AlleleProbs
    for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
      SampleAlleleFreqs(i, 1 );
      Loci(i)->SetAlleleProbs(Freqs[i]);
      if( Loci(i)->GetNumberOfLoci() > 1 )
	Loci(i)->SetHapPairProbs();
    }

    // Sample for allele frequency dispersion parameters, eta, using
    // Metropolis random-walk.
    if(  IsHistoricAlleleFreq ){
      Number++;
      for( int k = 0; k < Populations; k++ ){
         double mineta = 0;
	 vector< Vector_d > munew;
	 // Sample eta from truncated log-normal distribution.
	do{
           etanew = exp( gennor( log( eta[k] ), etastep[k] ) );
	}while( etanew > 5000.0 );
	// Prior log-odds ratio         
	LogPostRatio = ( psi[k] - 1 ) * (log(etanew) - log(eta[k]))
           - tau[k] * ( etanew - eta[k] );
	// Log-likelihood ratio; numerator of integrating constant
	LogPostRatio += 2 * Loci.GetNumberOfCompositeLoci()
           * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( eta[k] ) );
	for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){

	  Vector_d mu = GetPriorAlleleFreqs(j,k);
	  //mineta is a lower bound for proposal etanew
           if( mineta < 0.1 * eta[k] / mu.MinimumElement() )
              mineta = 0.1 * eta[k] / mu.MinimumElement();

	   //rescale munew so munew sums to etanew
           munew.push_back( mu * etanew / eta[k] );
           //Vector_d SumLogFreqs = GetStatsForEta(j, k );
	   double *SumLogFreqs = GetStatsForEta(j,k);
           for( int l = 0; l < Loci(j)->GetNumberOfStates(); l++ ){
	    // Denominator of integrating constant
              LogPostRatio += 2*(gsl_sf_lngamma( mu(l) ) - gsl_sf_lngamma( munew[j](l) ));
	    // SumLogFreqs = log phi_1 + log phi_2
              LogPostRatio += (munew[j](l) - mu(l))*SumLogFreqs[l];
           }
	   delete[] SumLogFreqs;
	}
	
	// Log acceptance probability = Log posterior ratio since the
	// proposal ratio (log-normal) cancels with prior.
	
	// Acceptance test.
	if( log( myrand() ) < LogPostRatio && mineta < etanew ){
           eta[k] = etanew;
           //UpdatePriorAlleleFreqs( k, munew );
           SumAcceptanceProb[k]++;
           NumberAccepted[k]++;
	}
	
	if( !( Number % w ) ){
	  etastep[k] = TuneEtaSampler[k].UpdateSigma( NumberAccepted[k] );
	  NumberAccepted[k] = 0;
	}
      }
      
      if( !( Number % w ) ){
	Number = 0;
      }
      
      if( iteration > BurnIn )
	for(int i=0;i<Populations;++i)SumEta[i]+=eta[i];
    }
    
    if( iteration > BurnIn && IsHistoricAlleleFreq ){
      UpdateFst();
    }
  }
}

/*
  Given a haplotype pair, h, and the ordered ancestry states at a locus,
  updates the counts of alleles observed in each state of ancestry.
 * should use hap pairs stored in Individual object
 */
void AlleleFreqs::UpdateAlleleCounts(int locus, int h[2], Vector_i ancestry )
{
  AlleleCounts[locus]( h[0], ancestry(1) )++;
  AlleleCounts[locus]( h[1], ancestry(0) )++;
}

//needs to be rewritten
void AlleleFreqs::UpdateAlleleCounts_HaploidData(int locus, std::vector<unsigned short >&genotype, int ancestry )
{
//     int xx;
//     if( Loci(locus)->GetNumberOfLoci() == 1 )
//       xx = genotype[0] - 1;
//     else{
//       Vector_i x = Loci(locus)->decodeGenotype(genotype);
//       xx = Loci(locus)->HapLoopGetDecimal( x );//fix this
//     }
//     AlleleCounts[locus]( xx, ancestry )++;

}

// get posterior mode of frequency of allele x, given locus and subpopulation
double AlleleFreqs::GetAlleleProbsMAP( int x, int ancestry , int locus)
{
   double P;
   if( x < Loci(locus)->GetNumberOfAllelesOfLocus(0) - 1 )
     // calculate posterior mode
     P = AlleleFreqsMAP[locus]( x, ancestry );
   else // frequency of last allele is set by subtracting sum of posterior modes of other alleles from 1 
   {
      P = 1;
      for( int j = 0; j < Loci(locus)->GetNumberOfAllelesOfLocus(0) - 1; j++ )
	P -= AlleleFreqsMAP[locus]( j, ancestry );
   }
   return P;
}


/* 
   gets probability of genotypes given ancestry states. 
   Probs is a KxK array in which rows and cols index 
   paternal and maternal locus ancestry states for each population.  
   this function will be redundant if GetGenotypeProbs is fixed to work with a haploid locus. can then call 
   GetGenotypeProbs directly
   Called by UpdateParameters method in Chromosome class and OnePopulationUpdate in Individual
*/
//TODO: haploid case
void AlleleFreqs::GetGenotypeProbs( double **Probs, int locus, std::vector<unsigned short >&genotype, 
				    std::vector<hapPair > &Haplotypes, bool diploid, bool fixed)
{
  if( diploid ){
    //here Probs should be a k x k array
    Loci(locus)->GetGenotypeProbs(Probs, Haplotypes, fixed, RandomAlleleFreqs);
  }
//   else{
//     // lines below should be replaced by a call to GetGenotypeProbs, which should be extended to 
//     // work with a haploid locus
//     //here Probs has a single column 
//     for(int i=0;i<Populations;++i)Probs[i][0]=0.0;
//      if( Loci(locus)->GetNumberOfLoci() == 1 ){
//         for( int pop = 0; pop < Populations; pop++ ){
// 	  Probs[pop][ 0 ] = GetAlleleProbs( genotype[0] - 1, pop , locus);
//         }
//      }
//      else{
//        Vector_i x = Loci(locus)->decodeGenotype(genotype);
//        int xx = Loci(locus)->HapLoopGetDecimal( x );
//         for( int pop = 0; pop < Populations; pop++ ){
// 	  Probs[pop][ 0 ] = GetAlleleProbs( xx - 1, pop , locus);
//         }
//      }
//   }
}

/**
 * Whether the object should remember the results of sampling.
 *
 * flag - integer representing a boolean. Set true (one) to remember
 *   sampled data. Set false (zero) during burn in.
 * i - locus at which to update
 */
void AlleleFreqs::SampleAlleleFreqs(int i, int flag )
{
  Vector_d freqs;
  
  for( int j = 0; j < Populations; j++ ){
    freqs = gendirichlet( PriorAlleleFreqs[i].GetColumn(j)
			  + AlleleCounts[i].GetColumn(j) );
    freqs.RemoveElement( Loci(i)->GetNumberOfStates() - 1 );
    Freqs[i].SetColumn( j, freqs );
    if( IsHistoricAlleleFreq ){
      freqs = gendirichlet( PriorAlleleFreqs[i].GetColumn(j)
			    + HistoricLikelihoodAlleleFreqs[i].GetColumn(j) );
      freqs.RemoveElement( Loci(i)->GetNumberOfStates() - 1 );
      HistoricAlleleFreqs[i].SetColumn( j, freqs );
    }
  }
  
  if( flag > 0 ){
    SumAlleleFreqs[i] += Freqs[i];
  }
}


/**
 * Sets the sum of allele frequencies over all sampling iterations to
 * zero.
 */
void AlleleFreqs::ResetSumAlleleFreqs()
{
for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
  SumAlleleFreqs[i].SetElements(0);
 }
}

void AlleleFreqs::SetMergedHaplotypes(Vector_d *alpha0, std::ofstream *LogFileStreamPtr, bool IsPedFile){
  //Note: alpha0 = alpha[0] in Latent
  for( int j = 0; j < GetNumberOfCompositeLoci(); j++ ){
    if( Loci(j)->GetNumberOfLoci() > 1 ){

      for( int k = 0; k < Populations; k++ )
	pp[k] = (*alpha0)(k) / alpha0->Sum();
      Loci(j)->SetDefaultMergeHaplotypes( pp, Freqs[j] );
      if(IsPedFile)
	*LogFileStreamPtr << "\"" << Loci(j)->GetLabel(0) << "\"" << endl;
      else
	*LogFileStreamPtr << Loci(j)->GetLabel(0) << endl;
      for( int k = 0; k < Loci(j)->GetNumberOfStates(); k++ ){
	*LogFileStreamPtr << k << " " << Loci(j)->GetMergedHaplotype(k) << endl;
      }
    }
  }
}

/**
 * Used when model is specified to allow allele freqs in admixed
 * population to vary from the "historical" allele freqs in the
 * unadmixed population.  Dirichlet distribution for allele freqs at
 * locus with k alleles is specified with (k - 1) frequency parameters
 * (mu) and with a single dispersion parameter (eta) this method
 * samples mu and eta, and updates Dirichlet parameters for the allele
 * frequencies
 */
double *AlleleFreqs::GetStatsForEta( int locus, int population)
{
  double *stats = new double[ Loci(locus)->GetNumberOfStates() ];
   for( int i = 0; i < Loci(locus)->GetNumberOfStates() - 1; i++ )
     stats[ i ] = log( Freqs[locus]( i, population ) ) + log( HistoricAlleleFreqs[locus]( i, population ) );
   stats[ Loci(locus)->GetNumberOfStates() - 1 ] = log( 1 - Freqs[locus].GetColumn( population ).Sum() )
     + log( 1 - HistoricAlleleFreqs[locus].GetColumn( population ).Sum() );
   return stats;
}

// void AlleleFreqs::UpdatePriorAlleleFreqs(int j, const vector<Vector_d>& mu)
// {
//   for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
//     PriorAlleleFreqs[i].SetColumn( j, mu[i] );
//     //    double sum;
//     //    Vector_d freqs;
//     //    for( int j = 0; j < Populations; j++ ){
//     //       freqs = PriorAlleleFreqs[i].GetColumn(j);
//     //       sum = freqs.Sum();
//     //       PriorAlleleFreqs[i].SetColumn( j, freqs * eta[j] / sum );
//     //    }
//   }
// }

void AlleleFreqs::SamplePriorAlleleFreqsMultiDim( int locus)
  // for a multi-allelic locus, we sample the Dirichlet prior allele freqs
  // by a Metropolis random walk
{
   vector<int> accept(Populations,0);
   Vector_d mu1, mu2;
   for( int j = 0; j < Populations; j++ ){
      double Proposal1=0, Proposal2=0, f1=0, f2=0;
      mu1 = PriorAlleleFreqs[locus].GetColumn(j) / eta[j];
      mu2 = gendirichlet( mu1 / MuProposal[locus][j].GetSigma() );
        
      for( int i = 0; i < Loci(locus)->GetNumberOfStates(); i++ ){
         f1 += 0.1 * log( mu1(i) ) + 0.1 * log( 1 - mu1(i) );
         f2 += 0.1 * log( mu2(i) ) + 0.1 * log( 1 - mu2(i) );
         Proposal1 += (eta[j] * mu2(i) - 1) * log( mu1(i) ) - gsl_sf_lngamma( eta[j] * mu2(i) );
         Proposal2 += (eta[j] * mu1(i) - 1) * log( mu2(i) ) - gsl_sf_lngamma( eta[j] * mu1(i) );
      }
      
      int numberofstates = mu1.GetNumberOfElements();
      for( int k = 0; k < numberofstates; k++ ){
         f1 -= Populations * gsl_sf_lngamma( mu1( k )* eta[j] );
         f2 -= Populations * gsl_sf_lngamma( mu2( k )* eta[j] );
	 f1 += gsl_sf_lngamma( mu1( k )* eta[j] + AlleleCounts[locus](k,j) );
	 f2 += gsl_sf_lngamma( mu2( k )* eta[j] + AlleleCounts[locus](k,j) );
	 f1 += gsl_sf_lngamma( mu1( k )* eta[j] + HistoricLikelihoodAlleleFreqs[locus](k,j) );
	 f2 += gsl_sf_lngamma( mu2( k )* eta[j] + HistoricLikelihoodAlleleFreqs[locus](k,j) );
      }
      if( log(myrand()) < f2 - f1 - Proposal2 + Proposal1 ){
	PriorAlleleFreqs[locus].SetColumn( j, mu2 * eta[j] );
         accept[j] = 1;
         MuProposal[locus][j].Event(true);
      }
      else
	MuProposal[locus][j].Event(false);
   }
}

void AlleleFreqs::SamplePriorAlleleFreqs1D( int locus)
  // for a diallelic locus, we can sample the proportion parameter of the beta distribution
  // conditional on the observed counts, the historic allele freqs and the dispersion parameter
  // using an adaptive rejection sampler 
{
   double lefttruncation = 0.1;
   double MuParameters[2];
   MatrixArray_i counts0(1);
   MatrixArray_d counts1(1);

// Construct adaptive rejection sampler for mu.
   counts0(0) = AlleleCounts[locus];
   counts1(0) = HistoricLikelihoodAlleleFreqs[locus];
   //warning message - move to wherever loci data are read
//    for(int row=0;row<counts0(0).GetNumberOfRows();++row)for(int col=0;col<counts0(0).GetNumberOfCols();++col)
//      if(counts0(0)(row,col)==0 && counts1(0)(row,col)==0){
//        cout<<"Warning: zero copies of allele ("<<row<<","<<col<<") in both admixed and unadmixed samples"<<endl;
//    //poss return name of comp locus
//    }
   DARS SampleMu( 0, 0, 0, MuParameters, 2, fMu, dfMu, ddfMu, counts0, counts1 );

   SampleMu.SetLeftTruncation( lefttruncation );
   for( int j = 0; j < Populations; j++ ){
//      sum = 0.0;
      MuParameters[0] = eta[j];
      MuParameters[1] = j;
//       MuParameters(1) = (float)AlleleCounts( NumberOfStates - 1, j );
//       MuParameters(2) = (float)HistoricLikelihoodAlleleFreqs( NumberOfStates - 1, j );
//       MuParameters(3) = (float)AlleleCounts( k, j );
//       MuParameters(4) = (float)HistoricLikelihoodAlleleFreqs( k, j );
      SampleMu.SetRightTruncation( eta[j] - lefttruncation );
      SampleMu.UpdateParameters( MuParameters, 2 );
      PriorAlleleFreqs[locus]( 0, j ) = SampleMu.Sample();
// Last prior frequency parameter is determined by; sum of mu's = eta.
      PriorAlleleFreqs[locus]( 1, j ) = eta[j] - PriorAlleleFreqs[locus]( 0, j );
   }
}


void AlleleFreqs::InitializeEtaOutputFile(AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log)
{
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
    if( options->getTextIndicator()  && options->getAnalysisTypeIndicator() >= 0)
      {
	//Dispersion parameters (eta)
	if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	  for( int k = 0; k < Populations; k++ ){
	    outputstream << "\"eta." << PopulationLabels[k].substr(1);
	  }
	}
	outputstream << endl;
      }
  }
}

void AlleleFreqs::OutputErgodicAvg( int samples,AdmixOptions *options, std::ofstream *avgstream)
{
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    for( int j = 0; j < Populations; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumEta[j] / samples << " ";
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
	  (*LogFileStreamPtr) << setprecision(6) << eta[j] << " ";
	}
      }
  //output to screen
    if( options->useCOUT() )
      {
	for( int j = 0; j < Populations; j++ ){
	  cout.width(9);
	  cout << setprecision(6) << eta[j] << " ";
	}
      }
  //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      for( int j = 0; j < Populations; j++ ){
	outputstream.width(9);
	outputstream << setprecision(6) << eta[j];
      }
     outputstream << endl;
    }
  }
}

Genome *AlleleFreqs::getLoci(){
  return &Loci;
}
CompositeLocus *AlleleFreqs::getLocus(int i){
  return (CompositeLocus *)(Loci(i));
}
int AlleleFreqs::GetNumberOfCompositeLoci(){
  return Loci.GetNumberOfCompositeLoci();
}
void AlleleFreqs::OutputAlleleFreqs(){
  if( IsRandom() ){
    allelefreqoutput->OutputAlleleFreqs(this);
  }
}
void AlleleFreqs::OutputFST(bool IsPedFile){
  for( int j = 0; j < GetNumberOfCompositeLoci(); j++ ){
    if(IsPedFile)
      fstoutputstream << "\"" << Loci(j)->GetLabel(0) << "\"";
    else
      fstoutputstream << Loci(j)->GetLabel(0);
    for(int k=0; k<Populations; ++k)fstoutputstream << " " << Fst[j][k];
    fstoutputstream << endl;
  }
}

void AlleleFreqs::Reset(){
  for( int i = 0; i < GetNumberOfCompositeLoci(); i++ ){
    Loci(i)->ResetScoreForMisSpecOfAlleleFreqs();
    AlleleCounts[i].SetElements(0);
  }
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


// should set allele freqs to the posterior mode
// but apparently just sets to current value
// probably doesn't matter where there is strong prior on allele freqs, and used in Chib algorithm 
// which just requires a value near the posterior mode
// set AlleleFreqsMAP and getAlleleFreqs are called by Individual object
void AlleleFreqs::setAlleleFreqsMAP()
{
  for(int i=0;i<GetNumberOfCompositeLoci();++i)
    AlleleFreqsMAP[i] = Freqs[i];
  //Loci(i)->setAlleleFreqsMAP(Freqs[i]);
}

/**
 * Indicates whether allele frequencies are fixed or random.
 *
 * Returns:
 * an integer representing a boolean. true (one) if allele frequencies
 * are random. false (zero) otherwise.
 * Should this return a boolean?
 */
int AlleleFreqs::IsRandom()
{
  return( RandomAlleleFreqs );
}

void AlleleFreqs::getLociCorrSummary(double *f[]){
  for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci();++i)
    f[0][i] = f[1][i] = LociCorrSummary[i];
}
/**PriorAlleleFreqs
 * Returns Dirichlet parameters for allele frequencies for a particular population and locus.
 * 
 * population - the number of the population (zero based)
 * locus - the number of the locus
 * returns:
 * a vector containing Dirichlet parameters for frequencies of each allele at the locus. 
 * Expected frequencies are calculated by dividing each parameter by the sum of parameters
 */
Vector_d AlleleFreqs::GetPriorAlleleFreqs( int locus, int population )
{
  return( PriorAlleleFreqs[locus].GetColumn( population ) );
}
/**
 * AlleleCounts is a matrix of counts of each allele in each
 * population this is the sufficient statistic for updating allele
 * frequencies.  Zero's the allele frequencies before they are updated.
 */
Matrix_i &AlleleFreqs::GetAlleleCounts(int locus)
{
  return( AlleleCounts[locus] );
}
Vector_i AlleleFreqs::GetAlleleCounts( int locus, int population )
{
  return( AlleleCounts[locus].GetColumn( population ) );
}
Vector_d AlleleFreqs::getAlleleFreqsMAP( int locus, int population )
{
  return( AlleleFreqsMAP[locus].GetColumn( population ) );
}
/**AlleleFreqs
 * Gets the frequencies of each haplotype in a composite locus.
 *
 * returns:
 * a matrix containing the frequencies of each allele
 * one row per allele state, one column per population
 */
Matrix_d &AlleleFreqs::GetAlleleFreqs(int locus)
{
  return( Freqs[locus] );
}
Matrix_d AlleleFreqs::GetSumAlleleFreqs(int locus)
{
  return( SumAlleleFreqs[locus] );
}

void AlleleFreqs::UpdateFst()
{
  for(int locus = 0; locus < GetNumberOfCompositeLoci(); ++locus){
    double q_admix,q_parental,f,H_admix, H_parental, H_combined, pbar;
    for( int k = 0; k < Populations; k++ ){
      H_admix = 0;
      H_parental = 0;
      H_combined = 0;
      
      for( int i = 0; i < Loci(locus)->GetNumberOfStates() - 1; i++ ){
	H_admix += Freqs[locus]( i, k ) * Freqs[locus]( i, k );
	H_parental += HistoricAlleleFreqs[locus]( i, k ) * HistoricAlleleFreqs[locus]( i, k );
	pbar = 0.5 * ( Freqs[locus]( i, k ) + HistoricAlleleFreqs[locus]( i, k ) );
	H_combined += pbar * pbar;
      }
      
      q_admix = 1 - Freqs[locus].GetColumn(k).Sum();
      H_admix += q_admix * q_admix;
      q_parental = 1 - HistoricAlleleFreqs[locus].GetColumn(k).Sum();
      H_parental += q_parental * q_parental;
      pbar = 0.5 * ( q_admix + q_parental );
      H_combined += pbar * pbar;
      
      H_combined = 1 - H_combined;
      H_admix = 1 - H_admix;
      H_parental = 1 - H_parental;
      f = ( H_combined - 0.5 * ( H_admix + H_parental ) ) / H_combined;
      Fst[locus][k] = 2*f / ( 1 + f );
    }
    for(int i=0;i < Populations;++i)
      SumFst[locus][i] += Fst[locus][i];
  }
}

double
fMu( Vector_d &parameters, MatrixArray_i& counts0, MatrixArray_d& counts1, double mu )
{
//    Vector_d rn(2), ri(2);
//    rn(0) = parameters(1);
//    rn(1) = parameters(2);
//    ri(0) = parameters(3);
//    ri(1) = parameters(4);
   int pop = (int)parameters(1);
   double eta = parameters(0);
   double prior = 0.1 * log( mu / eta ) + 0.1 * log( 1 - mu / eta );
   double f = prior - 2 * gsl_sf_lngamma( mu ) - 2 * gsl_sf_lngamma( eta - mu );
//   for( int i = 0; i < 2; i++ )
//      f += gsl_sf_lngamma( mu + ri(i) ) + gsl_sf_lngamma( eta - mu + rn(i) );
   f += gsl_sf_lngamma( mu+counts0(0)(0,pop) ) + gsl_sf_lngamma( eta-mu+counts0(0)(1,pop) );
   f += gsl_sf_lngamma( mu+counts1(0)(0,pop) ) + gsl_sf_lngamma( eta-mu+counts1(0)(1,pop) );

   return f;
}

double
dfMu( Vector_d &parameters, MatrixArray_i& counts0, MatrixArray_d& counts1, double mu )
{
//    Vector_d rn(2), ri(2);
//    rn(0) = parameters(1);
//    rn(1) = parameters(2);
//    ri(0) = parameters(3);
//    ri(1) = parameters(4);
   int pop = (int)parameters(1);
   double eta = parameters(0), x, y1, y2;
   double prior = 0.1 / mu - 0.1 / ( eta - mu );
   double f = prior;
   x = parameters(0) - mu;
   if(mu < 0)cout<<"\nError in dfMu in compositelocus.cc - arg mu to ddigam is negative\n"; 
   ddigam( &mu, &y1 );
   if(x < 0)cout<<"\nError in dfMu in compositelocus.cc - arg x to ddigam is negative\n"; 
   ddigam( &x, &y2 );
   f += 2 * ( y2 - y1 );

//    for( int i = 0; i < 2; i++ ){
//       x = mu + ri(i);
//       ddigam( &x, &y2 );
//       f += y2;
//       x = eta - mu + rn(i);
//       ddigam( &x, &y2 );
//       f -= y2;
//    }

   x = mu + counts0(0)(0,pop);
   ddigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts0(0)(1,pop);
   ddigam( &x, &y2 );
   f -= y2;

   x = mu + counts1(0)(0,pop);
   ddigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts1(0)(1,pop);
   ddigam( &x, &y2 );
   f -= y2;

   return f;
}

double
ddfMu( Vector_d &parameters, MatrixArray_i& counts0, MatrixArray_d& counts1, double mu )
{
//    Vector_d rn(2), ri(2);
//    rn(0) = parameters(1);
//    rn(1) = parameters(2);
//    ri(0) = parameters(3);
//    ri(1) = parameters(4);
   int pop = (int)parameters(1);
   double eta = parameters(0), x, y1, y2;
   double prior = -0.1 / (mu*mu) - 0.1 / (( eta - mu ) * ( eta - mu ) );
   double f = prior;
   x = parameters(0) - mu;
   trigam( &mu, &y1 );
   trigam( &x, &y2 );
   f -= 2 * ( y2 + y1 );

//    for( int i = 0; i < 2; i++ ){
//       x = mu + ri(i);
//       trigam( &x, &y2 );
//       f += y2;
//       x = eta - mu + rn(i);
//       trigam( &x, &y2 );
//       f += y2;
//    }

   x = mu + counts0(0)(0,pop);
   trigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts0(0)(1,pop);
   trigam( &x, &y2 );
   f += y2;

   x = mu + counts1(0)(0,pop);
   trigam( &x, &y2 );
   f += y2;
   x = eta - mu + counts1(0)(1,pop);
   trigam( &x, &y2 );
   f += y2;

   return f;
}

