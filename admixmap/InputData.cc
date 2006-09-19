/**
 *   ADMIXMAP
 *   InputData.cc 
 *   Class to read and check all input data files
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "InputData.h"
#include "AdmixOptions.h"
#include "utils/StringConvertor.h"
#include "utils/DataReader.h"
#include "Genome.h"
//#include "Chromosome.h"
#include "Comms.h"

#include <string>
#include <sstream>

using namespace std;

///Extracts population labels from header line of allelefreq input file
void InputData::getPopLabels(const Vector_s& data, size_t Populations, Vector_s& labels)
{
  if(data.size() != Populations+1){cout << "Error in getPopLabels\n";exit(1);}

  for (size_t i = 0; i < Populations; ++i) {
    labels.push_back( StringConvertor::dequote(data[i+1]) );
  }
}
void getLabels(const Vector_s& data, string *labels)
{
  for (size_t i = 0, index = 0; i < data.size(); ++i) {
    labels[index++] = StringConvertor::dequote(data[i]);
  }
}

#ifdef PARALLEL
#include "utils/StringSplitter.h"
void InputData::readGenotypesFile(const char *fname, Matrix_s& data)
{
  int worker_rank = MPI::COMM_WORLD.Get_rank() - 2;
  int NumWorkers = MPI::COMM_WORLD.Get_size() - 2;

  if (0 == fname || 0 == strlen(fname)) return;

  ifstream in(fname);
  if (!in.is_open()) {
    string msg = "Cannot open file for reading: \"";
    msg += fname;
    msg += "\"";
    throw runtime_error(msg);
  }

  data.clear();
  try {
    StringSplitter splitter;
    string line;        
    int linenumber = 0;
    Vector_s empty;

    while (getline(in, line)) {

      if( ( ( (linenumber-1)%NumWorkers) == worker_rank) || linenumber==0){
	if (!StringConvertor::isWhiteLine(line.c_str())) {//skip blank lines
	  data.push_back(splitter.split(line.c_str()));//split lines into strings
	}
      }
      else data.push_back(empty);//insert empty string vector

      ++linenumber;
    }
  } catch (...) {
    in.close();
    throw;
  }
}
#endif

InputData::InputData()
{
}

InputData::~InputData()
{
}

void InputData::readData(AdmixOptions *options, LogWriter &Log)
{
  Log.setDisplayMode(Quiet);
  try
    {
      // Read all input files.
      DataReader::ReadData(options->getLocusFilename(), locusData_, Log);   //locusfile
      DataReader::convertMatrix(locusData_, locusMatrix_, 1, 1,2);//drop first row, first col and last col

#ifdef PARALLEL
      if(Comms::isMaster()) Log << "Loading " << options->getGenotypesFilename() << ".\n";
      if(Comms::isWorker())//only workers read genotypes
	readGenotypesFile(options->getGenotypesFilename(), geneticData_); //genotypes file
#else
      DataReader::ReadData(options->getGenotypesFilename(), geneticData_, Log); 
#endif
      if(Comms::isMaster() || Comms::isWorker()){
	DataReader::ReadData(options->getCovariatesFilename(), inputData_, covariatesMatrix_,Log);     //covariates file
	DataReader::ReadData(options->getOutcomeVarFilename(), outcomeVarData_,outcomeVarMatrix_, Log);//outcomevar file
	DataReader::ReadData(options->getCoxOutcomeVarFilename(), coxOutcomeVarData_, Log);            //coxoutcomevar file
	DataReader::convertMatrix(coxOutcomeVarData_, coxOutcomeVarMatrix_, 1, 0,0);//drop first row in conversion

      }
      if(Comms::isFreqSampler()){//only one process reads freq files
	DataReader::ReadData(options->getAlleleFreqFilename(), alleleFreqData_, Log);
	DataReader::ReadData(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_, Log);            
	DataReader::ReadData(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_, Log);
	DataReader::ReadData(options->getEtaPriorFilename(), etaPriorData_,etaPriorMatrix_,  Log);
      }
      DataReader::ReadData(options->getReportedAncestryFilename(), reportedAncestryData_, reportedAncestryMatrix_, Log);

      Log << "\n";
      
    } catch (const exception& e) {
    cerr << "\nException occured during parsing of input file: \n" << e.what() << endl;
    exit(1);
  }
  catch(string s){
    cerr << "\nException occured during parsing of input file: \n" << s << endl;;
    exit(1);
  }
  NumSimpleLoci = getNumberOfSimpleLoci();
  NumCompositeLoci = determineNumberOfCompositeLoci();
  if(Comms::isWorker()) NumIndividuals = geneticData_.size() - 1;

  CheckData(options, Log);
}

void InputData::CheckData(AdmixOptions *options, LogWriter &Log){
  Log.setDisplayMode(Quiet);
  if(Comms::isWorker())
    {
      IsPedFile = determineIfPedFile();
      CheckGeneticData(options);
    }

  double threshold = 100.0;if(options->getHapMixModelIndicator())threshold /= options->getRhoPriorMean();
  checkLocusFile(options->getgenotypesSexColumn(), threshold, options->CheckData());
  //locusMatrix_ = locusMatrix_.SubMatrix(1, locusMatrix_.nRows() - 1, 1, 2);//remove header and first column of locus file

  if(Comms::isFreqSampler())CheckAlleleFreqs(options, Log);
#ifdef PARALLEL
  if(strlen(options->getAlleleFreqFilename()) || strlen(options->getPriorAlleleFreqFilename()) || strlen(options->getHistoricalAlleleFreqFilename())){
    //broadcast number of populations/block states, if inferred from file, from rank1
    int K = options->getPopulations();
    MPI::COMM_WORLD.Bcast(&K, 1, MPI::INT, 1);
    options->setPopulations(K);
  }

  int rank = MPI::COMM_WORLD.Get_rank();
  //first worker tells  nonworkers how many individuals there are
  if(rank==2){
    MPI::COMM_WORLD.Send(&NumIndividuals, 1, MPI::INT, 0, 0);//tell master
    MPI::COMM_WORLD.Send(&NumIndividuals, 1, MPI::INT, 1, 1);//tell freqsampler
  }
  else if(rank<2){
    MPI::Status status;
    MPI::COMM_WORLD.Recv(&NumIndividuals, 1, MPI::INT, 2, rank, status);
  }
#endif
  ReadPopulationLabels(options);

  if(Comms::isMaster() || Comms::isWorker() ){
    //detects regression model
    if(strlen( options->getOutcomeVarFilename() ) || strlen( options->getCoxOutcomeVarFilename() )){//if outcome specified
      if ( strlen( options->getOutcomeVarFilename() ) != 0 )
	CheckOutcomeVarFile( options, Log);
      if ( strlen( options->getCoxOutcomeVarFilename() ) != 0 ){
	OutcomeType.push_back( CoxData );
	if(options->CheckData())
	  CheckCoxOutcomeVarFile( Log);
      }
      if ( strlen( options->getCovariatesFilename() ) != 0 )
	CheckCovariatesFile(Log);
      //append population labels to covariate labels
      if(!options->getHapMixModelIndicator() && !options->getTestForAdmixtureAssociation()){
	for( vector<string>::const_iterator i = PopulationLabels.begin()+1; i !=PopulationLabels.end(); ++i ){
	  CovariateLabels.push_back("slope." + *i); 
	}
      }
    }
  
    if ( strlen( options->getReportedAncestryFilename() ) != 0 )
	CheckRepAncestryFile(options->getPopulations(), Log);
  }

  if(NumIndividuals > 1){
    Log.setDisplayMode(Quiet);
    Log << NumIndividuals << " individuals\n";
  }
}
///determine number of individuals by counting lines in genotypesfile 
int InputData::getNumberOfIndividuals()const {
  return(NumIndividuals);
}

///determine number of loci by counting rows of locusfile
int InputData::getNumberOfSimpleLoci()const {
  return(locusData_.size() - 1);
}
///determines number of composite loci from locusfile
unsigned InputData::determineNumberOfCompositeLoci()const{
  unsigned NumberOfCompositeLoci = locusMatrix_.nRows();
  for( unsigned i = 0; i < locusMatrix_.nRows(); i++ )
    if( !locusMatrix_.isMissing(i,1) && locusMatrix_.get( i, 1 ) == 0.0 ) NumberOfCompositeLoci--;
  return NumberOfCompositeLoci;
}
/// Determine if genotype table is in pedfile format by testing if number of strings in row 1 equals
/// twice the number of strings in the header row minus one. 
/// 
bool InputData::determineIfPedFile()const {
  const bool isPedFile = (bool)(2*geneticData_[0].size() - 1 == geneticData_[1].size());

  return (isPedFile);
}

///checks number of loci in genotypes file is the same as in locusfile, 
///determines if there is a sex column
/// and each line of genotypesfile has the same number of cols.
void InputData::CheckGeneticData(AdmixOptions *options)const{
  const size_t numLoci = locusData_.size() - 1; //number of loci in locus file
  int SexCol = 0;
  // Determine if "Sex" column present in genotypes file.
  if (numLoci == geneticData_[0].size() - 1) {
    SexCol = 0;//no sex col
  } else if (numLoci == geneticData_[0].size() - 2) {
    SexCol  = 1;//sex col
  } else {//too many cols
    cerr << "Error: " << numLoci << " loci in locus file but " <<  geneticData_[0].size() - 1 << " loci in genotypes file." << endl;
    exit(2);
  }
  options->setgenotypesSexColumn(SexCol);

  if(options->CheckData()){
    unsigned ExpCols;
#ifdef PARALLEL
    const int rank = MPI::COMM_WORLD.Get_rank()-2;
    const int numworkers = MPI::COMM_WORLD.Get_size()-2;
#else
    const int rank = 0;
    const int numworkers = 1;
#endif
    for(int i = rank + 1; i <= NumIndividuals; i+=numworkers){
      //should use logmsg
      if (IsPedFile) 
	ExpCols = 2*NumSimpleLoci + SexCol;
      else
	ExpCols = NumSimpleLoci + SexCol;
      
      if (geneticData_[i].size()-1 != ExpCols) {//check each row of genotypesfile has the right number of fields
	cerr << "Wrong number of entries ("<< geneticData_[i].size() <<")  in line "<<i+1<<" of genotypesfile" << endl;
	exit(1);
      }
    }
  }
}

bool InputData::distancesAreInCentiMorgans()const{
  bool distancesincM  = false;
  string distance_header = locusData_[0][2];
  if(distance_header.find("cm")!=string::npos || distance_header.find("CM")!=string::npos 
     || distance_header.find("cM")!=string::npos) 
    distancesincM = true;
  return distancesincM;
}

void InputData::checkLocusFile(int sexColumn, double threshold, bool check){
  // if check = true, Checks that loci labels in locusfile are unique and that they match the names in the genotypes file.
  //also extracts locus labels
  bool flag = false;

  if(distancesAreInCentiMorgans())threshold *= 100.0;
  for (size_t i = 1; i < locusData_.size(); ++i) {//rows of locusfile
    if(check && Comms::isFreqSampler()){
      //check number of alleles is >1
      if(locusMatrix_.get(i-1,0) <2){
	cerr << "ERROR on line " << i+1 << " of locusfile: number of alleles must be >1." << endl;
	flag = true;
      }

      //check distances are not negative
      if(locusMatrix_.get(i-1,1) < 0.0){
	flag = true;
	cerr<<"Error: distance on line "<<i<<" of locusfile is negative."<<endl;
      }
      //check distances are not too large 
      if(locusMatrix_.get(i-1,1) >= threshold) {
	//flag = true;
	if(locusMatrix_.get(i-1,1) != 100 )//for backward-compatibility; no warning if 100 used to denote new chromosome      
	  cerr << "Warning: distance of " <<locusMatrix_.get(i-1,1)<< "  at locus " <<i<<endl;
	locusMatrix_.isMissing(i-1,1, true);//missing value for distance denotes new chromosome
      }
      
      // Check loci names are unique    
      for (size_t j = 0; j < i-1; ++j) {   
	if (locusData_[i][0] == locusData_[j][0]) {
	  flag = true;
	  cerr << "Error in locusfile. Two different loci have the same name: "
	       << locusData_[i][0] << endl;
	}
      }

    }//end if check
    LocusLabels.push_back(StringConvertor::dequote(locusData_[i][0]));
  }//end loop over loci
  if(flag)exit(1);
  if(check && Comms::isWorker()){
    const size_t numLoci = locusData_.size() - 1;//number of simple loci

    // Compare loci names in locus file and genotypes file.
    for (size_t i = 1; i <= numLoci; ++i) {
      if (StringConvertor::dequote(locusData_[i][0]) != StringConvertor::dequote(geneticData_[0][i + sexColumn])) {
	cout << "Error. Locus names in locus file and genotypes file are not the same." << endl;
	cout << "Locus names causing an error are: " << locusData_[i][0] << " and " 
	     << geneticData_[0][i + sexColumn] << endl;
	//cout << options->getgenotypesSexColumn() << endl;
	exit(2);
      }
    }
  } 
}
 
void InputData::ReadPopulationLabels(AdmixOptions *options){
  if(strlen(options->getAlleleFreqFilename()) || strlen(options->getPriorAlleleFreqFilename()) || strlen(options->getHistoricalAlleleFreqFilename())){
    if(strlen(options->getAlleleFreqFilename()))
      DataReader::ReadHeader(options->getAlleleFreqFilename(), PopulationLabels);
    else if(strlen(options->getPriorAlleleFreqFilename()))
      DataReader::ReadHeader(options->getPriorAlleleFreqFilename(), PopulationLabels);
    else if(strlen(options->getHistoricalAlleleFreqFilename()))
      DataReader::ReadHeader(options->getHistoricalAlleleFreqFilename(), PopulationLabels);

  }
  else{
    //set default pop labels
    for( int j = 0; j < options->getPopulations(); j++ ){
      stringstream poplabel;
      if(options->getHapMixModelIndicator()) poplabel << "BlockState" << j+1;
      else poplabel << "Pop" << j+1;
      PopulationLabels.push_back(poplabel.str());
    }
  }
}
 
////checks consistency of supplied allelefreqs with locusfile
///and determines number of populations and population labels.
void InputData::CheckAlleleFreqs(AdmixOptions *options, LogWriter &Log){
  string freqtype = "";
  bool infile = false;//indicates whether either of the three allelefreq files are specified
  int nrows=0, expectednrows=0;
  int Populations = options->getPopulations();
  int NumberOfStates = 0;

  unsigned index = 0;
  for(unsigned i = 0; i < NumCompositeLoci; ++i){
    int states = 1;
    do{
      states *= (int)locusMatrix_.get( index, 0 );
      index++;
    }
    while( index < locusMatrix_.nRows() - 1 && !locusMatrix_.isMissing(index, 1) && locusMatrix_.get( index, 1 ) == 0 );
    NumberOfStates += states;
  }


  //fixed allele freqs
  if( strlen( options->getAlleleFreqFilename() ) ){
    freqtype = "";
    infile = true;
    nrows = alleleFreqData_.size()-1;
    expectednrows = NumberOfStates-NumCompositeLoci;
    Populations = alleleFreqData_[0].size() - 1;// -1 for ids in first col
    //getPopLabels(alleleFreqData_[0], Populations, PopulationLabels);
  }
  
  //Historic allelefreqs
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    freqtype = "historic";
    infile = true;
    nrows = historicalAlleleFreqData_.size();
    expectednrows = NumberOfStates+1;
    Populations = historicalAlleleFreqData_[0].size() - 1;
    //getPopLabels(historicalAlleleFreqData_[0], Populations, PopulationLabels);

  }
  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
    freqtype = "prior";
    infile = true;
    nrows = priorAlleleFreqData_.size();
    expectednrows = NumberOfStates+1;
    Populations = priorAlleleFreqData_[0].size() - 1;
    //getPopLabels(priorAlleleFreqData_[0], Populations, PopulationLabels);
  }
  if(infile){
    if(nrows != expectednrows){
      Log << "Incorrect number of rows in " << freqtype << "allelefreqfile.\n" 
	  << "Expecting " << expectednrows << " rows, but there are " << nrows << " rows.\n";
      exit(1);
    }
    options->setPopulations(Populations);
  }
  else{//'populations' option
    if(Populations < 1){
      Log << "ERROR: populations = " << options->getPopulations() << "\n";
      exit(1);
    }

    //     for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    //       if(Loci->GetNumberOfStates(i) < 2){
    // 	Log << "ERROR: The number of alleles at a locus is " << Loci->GetNumberOfStates(i) << "\n";
    // 	exit(1);
    //       }
    //     }
  }
}

void InputData::CheckOutcomeVarFile(AdmixOptions* const options, LogWriter& Log){
  //check outcomevarfile and genotypes file have the same number of rows
  if( (int)outcomeVarMatrix_.nRows() - 1 != (NumIndividuals - options->getTestOneIndivIndicator()) ){
    stringstream s;
    s << "ERROR: Genotypes file has " << NumIndividuals << " observations and Outcomevar file has "
      << outcomeVarMatrix_.nRows() - 1 << " observations.\n";
    throw(s.str());
  }
  //check the number of outcomes specified is not more than the number of cols in outcomevarfile
  int Firstcol = options->getTargetIndicator();
  int NumOutcomes = options->getNumberOfOutcomes();
  if(strlen(options->getCoxOutcomeVarFilename())) --NumOutcomes;

  int numoutcomes = NumOutcomes;
  if(NumOutcomes > -1){//options 'numberofregressions' used
    if((int)outcomeVarMatrix_.nCols() - Firstcol < NumOutcomes){
      numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;//adjusts if too large
      Log << "ERROR: 'outcomes' is too large, setting to " << numoutcomes;
    }
  }
  else numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;
  options->setNumberOfOutcomes(numoutcomes);

  RegressionType RegType = None;  
  if(numoutcomes >0){
    //RegType = Both;
    //extract portion of outcomevarfile needed
    std::string* OutcomeVarLabels = new string[ outcomeVarMatrix_.nCols() ];
    getLabels(outcomeVarData_[0], OutcomeVarLabels);
    DataMatrix Temp = outcomeVarMatrix_.SubMatrix(1, NumIndividuals, Firstcol, Firstcol+numoutcomes-1);
    outcomeVarMatrix_ = Temp;
    
    //determine type of outcome - binary/continuous
    for( int j = 0; j < numoutcomes; j++ ){
      OutcomeType.push_back( Binary );
      for(int i = 0; i < NumIndividuals; ++i)
	if(!outcomeVarMatrix_.isMissing(i, j) && !(outcomeVarMatrix_.get( i, j ) == 0 || outcomeVarMatrix_.get( i, j ) == 1) ){
	  OutcomeType[j] =  Continuous ;
	  break;
	}
      //in this way, the outcome type is set as binary only if all individuals have outcome values of 1 or 0
      //otherwise, a continuous outcome of 1.0 or 0.0 could lead to the type being wrongly set to binary.
      
      //need to check for allmissing
      //     if(i == NumIndividuals){
      //       Log << "ERROR: all outcomes missing\n";
      //       exit(1);
      //     }
      
      Log << "Regressing on ";    
      if( OutcomeType[j] == Binary ){
	Log << "Binary variable: ";
	if(numoutcomes==1)RegType = Logistic;//one logistic regression
	else {
	  if(RegType == Logistic) RegType = Mlogistic;//more than one logistic regression
	  else RegType = Multiple;//linear and logistic 
	}
      }
      else if(OutcomeType[j] == Continuous ){
	Log << "Continuous variable: ";
	if(numoutcomes==1)RegType = Linear;//one linear regression
	else {
	  if(RegType == Linear) RegType = Mlinear;//more than one linear regression
	  else RegType = Multiple;//linear and logistic 
	}
      }
      Log << outcomeVarData_[0][j+Firstcol];
      Log << ".\n";
      OutcomeLabels.push_back(outcomeVarData_[0][j+Firstcol]);
    }
    Log << "\n";
    delete[] OutcomeVarLabels;
  }
  options->setRegType(RegType);
}

void InputData::CheckCoxOutcomeVarFile(LogWriter &Log)const{
  if(coxOutcomeVarMatrix_.nCols() !=3){
    Log << "ERROR: 'coxoutcomevarfile should have 3 columns but has " << coxOutcomeVarMatrix_.nCols() << "\n";
    exit(1);
  }
  if( (int)coxOutcomeVarMatrix_.nRows() != NumIndividuals ){
    stringstream s;
    s << "ERROR: Genotypes file has " << NumIndividuals << " observations and coxoutcomevar file has "
      << coxOutcomeVarMatrix_.nRows() - 1 << " observations.\n";
    throw(s.str());
  }

  bool flag = false;
  for(unsigned i = 0; i < coxOutcomeVarMatrix_.nRows(); ++i)
    if(coxOutcomeVarMatrix_.get(i, 0) >= coxOutcomeVarMatrix_.get(i, 1)){
      Log << "Error in coxoutcomevarfile on line " <<i << " : finish times must be later than start times\n";
      flag = true;
    }
  if(flag)exit(1);

}

void InputData::CheckCovariatesFile(LogWriter &Log){
  if( NumIndividuals != (int)covariatesMatrix_.nRows() - 1 ){
    Log << "ERROR: Genotypes file has " << NumIndividuals << " observations and Covariates file has "
	<< covariatesMatrix_.nRows() - 1 << " observations.\n";
    exit(1);
  }
  for (size_t i = 0; i < inputData_[0].size(); ++i) {
    CovariateLabels.push_back(StringConvertor::dequote(inputData_[0][i]));
  }
}

void InputData::CheckRepAncestryFile(int populations, LogWriter &Log)const{
  if( (int)reportedAncestryMatrix_.nRows() != 2 * NumIndividuals ){
    Log << "ERROR: " << "ReportedAncestry file has " << reportedAncestryMatrix_.nRows() << " rows\n"
	<<"Genotypesfile has " << NumIndividuals << " rows\n";
    exit(1);}
  if( (int)reportedAncestryMatrix_.nCols() != populations ){
    Log << "ERROR: " << "ReportedAncestry file has " << reportedAncestryMatrix_.nCols() << " cols\n"
	<< "AlleleFreq file has "<< populations << " cols\n";
    exit(1);
  }
}

///determines if an individual is female
bool InputData::isFemale(int i)const{
  //if (options->getgenotypesSexColumn() == 1) {
  int sex = StringConvertor::toInt(geneticData_[i][1]);
  if (sex > 2) {
    cout << "Error: sex must be coded as 0 - missing, 1 - male or 2 - female.\n";
    exit(0);
  }        
  //}
  return (bool)(sex==2);
}

vector<unsigned short> InputData::GetGenotype(unsigned locus, int individual, int SexColumn)const{
  vector<unsigned short> g;
  bool isXlocus = false;
  if(locusData_[0].size()==4){
    const string s1("X"), s2("x");
    const string chrmlabel = StringConvertor::dequote(locusData_[locus+1][3]);
    isXlocus = ( (chrmlabel == s1) || (chrmlabel == s2) );
  }
  int col = 1 + SexColumn + locus;
  if (IsPedFile)col = 1 + SexColumn + 2*locus;
  //strip quotes from string
  const std::string str = StringConvertor::dequote(geneticData_[individual][col]);
  //look for , or / 
  string::size_type i = str.find_first_of(",/");
  //extract first allele as portion of string up to first ',' or '/'
  //NOTE: if string consists only of ',' or '/' or anything non-numeric, genotype is taken as missing
  g.push_back(atoi(str.substr(0,i).c_str()));

  if(!isXlocus || isFemale(individual)){//expect diploid genotype
    if( i != string::npos){// , or / found
      //extract second allele as portion of string after first ',' or '/'
      // NOTE: if nothing after, allele is taken as 0
      g.push_back(atoi(str.substr(i+1,str.length()-i).c_str()));
    }
    else{
      //if empty string, interpret as missing genotype for backward compatibility
      if(str.length()==0)g.push_back(0);
      else {//string with only one int
	cerr << "Error in genotypesfile: expected diploid genotype for Individual " << individual << " at locus " << locus 
	     << " but found haploid genotype"<< endl;
	exit(1);
      }
    }
  }
  //   else{
  //   //check male X genotypes are haploid
  //     if(i != string::npos){//found another allele
  // 	cerr << "Error in genotypesfile: expected haploid genotype for Individual " << individual << " at locus " << locus 
  // 	     << " but found diploid genotype"<< endl;
  //     }
  //   }   
  return g;  
}

void InputData::GetGenotype(int i, int SexColumn, const Genome &Loci, vector<genotype>* genotypes, bool** Missing)const{
  unsigned int simplelocus = 0;//simple locus counter
  unsigned complocus = 0;
  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    for(unsigned int j = 0; j < Loci.GetSizeOfChromosome(c); ++j){
      genotype G;
      // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype
#ifdef PARALLEL
      const int numLoci = 1; 
#else
      const int numLoci = Loci(complocus)->GetNumberOfLoci();
#endif
      
      unsigned int count = 0;
      for (int locus = 0; locus < numLoci; locus++) {
#ifdef PARALLEL
	const int numalleles = 2;
#else
	const int numalleles = Loci(complocus)->GetNumberOfAllelesOfLocus(locus);
#endif
	vector<unsigned short> g = GetGenotype(simplelocus, i, SexColumn);
	if(g.size()==2)
	  if( (g[0] > numalleles) || (g[1] > numalleles))
	    throwGenotypeError(i, simplelocus, Loci(complocus)->GetLabel(0), 
			       g[0], g[1], numalleles );
	  else if (g.size()==1)
	    if( (g[0] > numalleles))
	      throwGenotypeError(i, simplelocus, Loci(complocus)->GetLabel(0), 
				 g[0], 0, numalleles );

	simplelocus++;
	G.push_back(g);
	count += g[0];
      }
      
      Missing[c][j] = (count == 0);
      
      genotypes->push_back(G);
      ++complocus;
    }
  }
}

void InputData::throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles)const{
  cerr << "Error in genotypes file:\n"
       << "Individual " << ind << " at locus " << locus <<" (" << label << ")"
       << " has genotype " << g0 << ", " << g1 << " \n"
       << "Number of allelic states at locus = " << numalleles << "\n";
  if(ind == NumIndividuals)
    exit(1);
}

void InputData::getOutcomeTypes(DataType* T)const{
  for(unsigned i = 0; i < OutcomeType.size(); ++i)
    T[i] = OutcomeType[i];
}
DataType InputData::getOutcomeType(unsigned i)const{
  return OutcomeType[i];
}
const Matrix_s& InputData::getLocusData() const
{
  return locusData_;
}

const Matrix_s& InputData::getGeneticData() const
{
  return geneticData_;
}

const Matrix_s& InputData::getInputData() const
{
  return inputData_;
}

const Matrix_s& InputData::getOutcomeVarData() const
{
  return outcomeVarData_;
}

const Matrix_s& InputData::getAlleleFreqData() const
{
  return alleleFreqData_;
}

const Matrix_s& InputData::getHistoricalAlleleFreqData() const
{
  return historicalAlleleFreqData_;
}

const Matrix_s& InputData::getPriorAlleleFreqData() const
{
  return priorAlleleFreqData_;
}

const Matrix_s& InputData::getEtaPriorData() const
{
  return etaPriorData_;
}

const Matrix_s& InputData::getReportedAncestryData() const
{
  return reportedAncestryData_;
}

const DataMatrix& InputData::getEtaPriorMatrix() const
{
  return etaPriorMatrix_;
}

const DataMatrix& InputData::getLocusMatrix() const
{
  return locusMatrix_;
}

// const DataMatrix& InputData::getAlleleFreqMatrix() const
// {
//     return alleleFreqMatrix_;
// }

// const DataMatrix& InputData::getHistoricalAlleleFreqMatrix() const
// {
//     return historicalAlleleFreqMatrix_;
// }

// const DataMatrix& InputData::getPriorAlleleFreqMatrix() const
// {
//     return priorAlleleFreqMatrix_;
// }

const DataMatrix& InputData::getOutcomeVarMatrix() const
{
  return outcomeVarMatrix_;
}
const DataMatrix& InputData::getCoxOutcomeVarMatrix() const
{
  return coxOutcomeVarMatrix_;
}

const DataMatrix& InputData::getReportedAncestryMatrix() const
{
  return reportedAncestryMatrix_;
}

const DataMatrix& InputData::getCovariatesMatrix() const
{
  return covariatesMatrix_;
}
const Vector_s& InputData::GetPopLabels() const{
  return PopulationLabels;
}
Vector_s InputData::getOutcomeLabels()const{
  return OutcomeLabels;
}
const Vector_s& InputData::getLocusLabels()const{
  return LocusLabels;
}
const Vector_s InputData::getCovariateLabels()const{
  return CovariateLabels;
}

void InputData::Delete(){
  //erase string matrices
  for(unsigned i = 0; i < locusData_.size(); ++i)
    locusData_[i].clear();
  locusData_.clear();
  for(unsigned i = 0; i < geneticData_.size(); ++i)
    geneticData_[i].clear();
  geneticData_.clear();
  for(unsigned i = 0; i < inputData_.size(); ++i)
    inputData_[i].clear();
  inputData_.clear();
  for(unsigned i = 0; i < outcomeVarData_.size(); ++i)
    outcomeVarData_[i].clear();
  outcomeVarData_.clear();
  for(unsigned i = 0; i < coxOutcomeVarData_.size(); ++i)
    coxOutcomeVarData_[i].clear();
  coxOutcomeVarData_.clear();
  for(unsigned i = 0; i < alleleFreqData_.size(); ++i)
    alleleFreqData_[i].clear();
  alleleFreqData_.clear();
  for(unsigned i = 0; i < priorAlleleFreqData_.size(); ++i)
    priorAlleleFreqData_[i].clear();
  priorAlleleFreqData_.clear();
  for(unsigned i = 0; i < historicalAlleleFreqData_.size(); ++i)
    historicalAlleleFreqData_[i].clear();
  historicalAlleleFreqData_.clear();
  for(unsigned i = 0; i < etaPriorData_.size(); ++i)
    etaPriorData_[i].clear();
  etaPriorData_.clear();
  for(unsigned i = 0; i < reportedAncestryData_.size(); ++i)
    reportedAncestryData_[i].clear();
  reportedAncestryData_.clear();

  //erase data matrices 
  locusMatrix_.clear();
  covariatesMatrix_.clear();
  outcomeVarMatrix_.clear();
  coxOutcomeVarMatrix_.clear();
  //alleleFreqMatrix_.clear();
  //historicalAlleleFreqMatrix_.clear();
  //priorAlleleFreqMatrix_.clear();
  etaPriorMatrix_.clear();
  reportedAncestryMatrix_.clear();
}
