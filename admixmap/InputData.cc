/**
 *   ADMIXMAP
 *   InputData.cc 
 *   Class to read and check all input data files
 *   Copyright (c) 2005, 2006 LSHTM
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
#include "InputData.h"
#include "AdmixOptions.h"
#include "StringSplitter.h"
#include "StringConvertor.h"
#include "Genome.h"
#include "Chromosome.h"

#include <string>
#include <sstream>

using namespace std;

//Extracts population labels from header line of allelefreq input file
void InputData::getPopLabels(const Vector_s& data, size_t Populations, string **labels)
{
  if(data.size() != Populations+1){cout << "Error in getPopLabels\n";exit(1);}
  *labels = new string[ Populations ];

    for (size_t i = 0; i < Populations; ++i) {
      (*labels)[i] = StringConvertor::dequote(data[i+1]);
    }
}
void getLabels(const Vector_s& data, string *labels)
{
  for (size_t i = 0, index = 0; i < data.size(); ++i) {
    labels[index++] = StringConvertor::dequote(data[i]);
  }
}

void InputData::readFile(const char *fname, Matrix_s& data, LogWriter &Log)
{
    if (0 == fname || 0 == strlen(fname)) return;

    ifstream in(fname);
    if (!in.is_open()) {
        string msg = "Cannot open file for reading: \"";
        msg += fname;
        msg += "\"";
        throw runtime_error(msg.c_str());
    }
    else {
      Log << "Loading " << fname << ".\n";
    }

    data.clear();
    try {
        StringSplitter splitter;

        string line;        

        while (getline(in, line)) {
	  if (!StringConvertor::isWhiteLine(line.c_str())) {
	      data.push_back(splitter.split(line.c_str()));
            }
        }
    } catch (...) {
        in.close();
        throw;
    }
}

/**
 *  Auxilary function that converts Matrix_s to DataMatrix
 */
static void convertMatrix(const Matrix_s& data, DataMatrix& m)
{       
    const size_t numRows = data.size();

    // If there are no rows, return empty matrix.
    if (0 == numRows) return;

    // Verify that all rows have same length.
    const size_t numCols = data[0].size();
    for (size_t i = 1; i < numRows; ++i) {
        if (numCols != data[i].size()) {
            throw runtime_error("Invalid row length");
        }
    }
    
    // Form matrix.
    m.setDimensions(numRows, numCols);
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            if (StringConvertor::isMissingValue(data[i][j])) {
	      m.isMissing(i, j, true);
            } else {
	      m.set(i, j, StringConvertor::toFloat(data[i][j]));
            }
        }
    }
}

/**
 *  InputData members.
 */

InputData::InputData()
{
  PopulationLabels = 0;
}

InputData::~InputData()
{
  delete[] PopulationLabels;
}

void InputData::readData(AdmixOptions *options, LogWriter &Log)
{
  Log.setDisplayMode(Quiet);
  try
    {
      // Read all input files.
      readFile(options->getLocusFilename(), locusData_, Log);   //locusfile
      readFile(options->getGenotypesFilename(), geneticData_, Log); //genotypes file
      readFile(options->getCovariatesFilename(), inputData_, Log);         //covariates file
      readFile(options->getOutcomeVarFilename(), outcomeVarData_, Log);       //outcomevar file
      readFile(options->getAlleleFreqFilename(), alleleFreqData_, Log);
      readFile(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_, Log);            
      readFile(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_, Log);
      readFile(options->getEtaPriorFilename(), etaPriorData_, Log);
      readFile(options->getReportedAncestryFilename(), reportedAncestryData_, Log);

      Log << "\n";
      // Form matrices.
      convertMatrix(locusData_, locusMatrix_);
      ::convertMatrix(outcomeVarData_, outcomeVarMatrix_);
      ::convertMatrix(inputData_,  covariatesMatrix_);
      ::convertMatrix(alleleFreqData_, alleleFreqMatrix_);
      ::convertMatrix(historicalAlleleFreqData_, historicalAlleleFreqMatrix_);
      ::convertMatrix(priorAlleleFreqData_, priorAlleleFreqMatrix_);
      ::convertMatrix(etaPriorData_, etaPriorMatrix_);
      ::convertMatrix(reportedAncestryData_, reportedAncestryMatrix_);
      
    } catch (const exception& e) {
    cerr << "Exception occured during parsing input file: \n" << e.what() << endl;
    exit(1);
  }
  NumSimpleLoci = getNumberOfSimpleLoci();
  NumCompositeLoci = determineNumberOfCompositeLoci();
  NumIndividuals = getNumberOfIndividuals();

  Log.setDisplayMode(Quiet);
  IsPedFile = determineIfPedFile( options );
  CheckGeneticData(options);

  double threshold = 100.0;if(options->getHapMixModelIndicator())threshold /= options->getRhoPriorMean();
  checkLocusFile(options->getgenotypesSexColumn(), threshold, options->CheckData());
  locusMatrix_ = locusMatrix_.SubMatrix(1, locusMatrix_.nRows() - 1, 1, 2);//remove header and first column of locus file
  if ( strlen( options->getOutcomeVarFilename() ) != 0 )
    options->setRegType( CheckOutcomeVarFile( options->getNumberOfOutcomes(), options->getTargetIndicator(), Log));
  if ( strlen( options->getCovariatesFilename() ) != 0 )
    CheckCovariatesFile(Log);//detects regression model
  if ( strlen( options->getReportedAncestryFilename() ) != 0 )
    CheckRepAncestryFile(options->getPopulations(), Log);
  CheckAlleleFreqs(options, Log);
 
  if(NumIndividuals > 1){
    Log.setDisplayMode(Quiet);
    Log << NumIndividuals << " individuals\n";
  }
}

//determine number of individuals by counting lines in genotypesfile 
int InputData::getNumberOfIndividuals()const {
  return(geneticData_.size() - 1);
}

//determine number of loci by counting rows of locusfile
int InputData::getNumberOfSimpleLoci()const {
  return(locusData_.size() - 1);
}
//determines number of composite loci from locusfile
unsigned InputData::determineNumberOfCompositeLoci()const{
  unsigned NumberOfCompositeLoci = locusMatrix_.nRows()-1;
    for( unsigned i = 1; i < locusMatrix_.nRows(); i++ )
      if( !locusMatrix_.isMissing(i,2) && locusMatrix_.get( i, 2 ) == 0.0 ) NumberOfCompositeLoci--;
    return NumberOfCompositeLoci;
}

bool InputData::determineIfPedFile(AdmixOptions *options)const {
  // Determine if genotype table is in pedfile format by testing if number of strings in row 1 equals
  // twice the number of strings in the header row minus one. 
  // 
  const bool isPedFile = (bool)(2*geneticData_[0].size() - 1 == geneticData_[1].size());
  options->IsPedFile(isPedFile);

  return (isPedFile);
}

//checks number of loci in genotypes file is the same as in locusfile
void InputData::CheckGeneticData(AdmixOptions *options)const{

  const size_t numLoci = locusData_.size() - 1; //number of loci in locus file
  int sexcol;
  // Determine if "Sex" column present in genotypes file.
  if (numLoci == geneticData_[0].size() - 1) {
    sexcol = 0;
  } else if (numLoci == geneticData_[0].size() - 2) {
    sexcol  = 1;
  } else {
    cerr << "Error. Number of loci in genotypes file does not match number in locus file." << endl;
    exit(2);
  }
  options->setgenotypesSexColumn(sexcol);

  unsigned ExpCols;
  for(int i = 1; i <= NumIndividuals; ++i){
    //should use logmsg
    if (IsPedFile) 
      ExpCols = 2*NumSimpleLoci + sexcol;
    else
      ExpCols = NumSimpleLoci + sexcol;
    
    if (geneticData_[i].size()-1 != ExpCols) {//check each row of genotypesfile has the right number of fields
      cerr << "Wrong number of entries ("<< geneticData_[i].size() <<")  in line "<<i+1<<" of genotypesfile" << endl;
      exit(1);
    }
    
  }
}

void InputData::checkLocusFile(int sexColumn, double threshold, bool check){
  // Check that loci labels in locusfile are unique and that they match the names in the genotypes file.
  //also determines number of chromosomes
  //NumChromosomes = 0;
  bool flag = false;
  for (size_t i = 1; i < locusData_.size(); ++i) {//rows of locusfile
      if(check){
    //check distances are not negative
    if(locusMatrix_.get(i,2) < 0.0){
      flag = true;
      cerr<<"Error: distance on line "<<i<<" of locusfile is negative."<<endl;
    }
    //check distances are not too large 
    if(locusMatrix_.get(i,2) >= threshold) {
      //flag = true;
      if(locusMatrix_.get(i,2) < 100 )//for backward-compatibility; no warning if 100 used to denote new chromosome      
	cerr << "Warning: distance of " <<locusMatrix_.get(i,2)<< "  at locus " <<i<<endl;
      locusMatrix_.isMissing(i,2, true);//missing value for distance denotes new chromosome
    }
      
    // Check loci names are unique    
    for (size_t j = 0; j < i-1; ++j) {   
      if (locusData_[i][0] == locusData_[j][0]) {
	flag = true;
	cerr << "Error in locusfile. Two different loci have the same name: "
	     << locusData_[i][0] << endl;
      }
    }
    //if(locusMatrix_.isMissing(i,2))++NumChromosomes;
      }
    LocusLabels.push_back(StringConvertor::dequote(locusData_[i][0]));
  }//end loop over loci
  if(flag)exit(1);
  if(check){
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

//checks consistency of supplied allelefreqs with locusfile
//and determines number of populations and population labels
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
    nrows = alleleFreqMatrix_.nRows()-1;
    expectednrows = NumberOfStates-NumCompositeLoci;
    Populations = alleleFreqMatrix_.nCols() - 1;// -1 for ids in first col
    getPopLabels(alleleFreqData_[0], Populations, &PopulationLabels);
  }
  
  //Historic allelefreqs
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    freqtype = "historic";
    infile = true;
    nrows = historicalAlleleFreqMatrix_.nRows();
    expectednrows = NumberOfStates+1;
    Populations = historicalAlleleFreqMatrix_.nCols() - 1;
    getPopLabels(historicalAlleleFreqData_[0], Populations, &PopulationLabels);

  }
  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
      freqtype = "prior";
      infile = true;
      nrows = priorAlleleFreqMatrix_.nRows();
      expectednrows = NumberOfStates+1;
      Populations = priorAlleleFreqMatrix_.nCols() - 1;
      getPopLabels(priorAlleleFreqData_[0], Populations, &PopulationLabels);
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
    PopulationLabels = new string[ Populations ];
    for( int j = 0; j < Populations; j++ ){
      stringstream poplabel;
      string result;
      if(options->getHapMixModelIndicator()) poplabel << "BlockState" << j+1;
      else poplabel << "Pop" << j+1;
      result = poplabel.str();
      PopulationLabels[j] = result;
    }
//     for( int i = 0; i < NumberOfCompositeLoci; i++ ){
//       if(Loci->GetNumberOfStates(i) < 2){
// 	Log << "ERROR: The number of alleles at a locus is " << Loci->GetNumberOfStates(i) << "\n";
// 	exit(1);
//       }
//     }
  }
}

RegressionType InputData::CheckOutcomeVarFile(int NumOutcomes, int Firstcol, LogWriter& Log){
  //check outcomevarfile and genotypes file have the same number of cols
  if( (int)outcomeVarMatrix_.nRows() - 1 != NumIndividuals ){
    Log << "ERROR: Genotypes file has " << NumIndividuals << " observations and Outcomevar file has "
	<< outcomeVarMatrix_.nRows() - 1 << " observations.\n";
    exit(1);
  }
  //check the number of outcomes specified is not more than the number of cols in outcomevarfile
  int numoutcomes = NumOutcomes;
  if(NumOutcomes > -1){//options 'numberofregressions' used
    if((int)outcomeVarMatrix_.nCols() - Firstcol < NumOutcomes){
      numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;//adjusts if too large
      Log << "ERROR: 'outcomes' is too large, setting to " << numoutcomes;
    }
  }
  else numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;

  RegressionType RegType = None;  
  if(numoutcomes >0){
    RegType = Both;
    //extract portion of outcomevarfile needed
    std::string* OutcomeVarLabels = new string[ outcomeVarMatrix_.nCols() ];
    getLabels(outcomeVarData_[0], OutcomeVarLabels);
    DataMatrix Temp = outcomeVarMatrix_.SubMatrix(1, NumIndividuals, Firstcol, Firstcol+numoutcomes-1);
    outcomeVarMatrix_ = Temp;
    
    //determine type of outcome - binary/continuous
    OutcomeType = new DataType[numoutcomes];
    for( int j = 0; j < numoutcomes; j++ ){
      OutcomeType[j] = Binary; // default type
      for(int i = 0; i < NumIndividuals; ++i)
	if(!outcomeVarMatrix_.isMissing(i, j) && 
	   !(outcomeVarMatrix_.get( i, j ) == 0 || outcomeVarMatrix_.get( i, j ) == 1) ) {
	  OutcomeType[j] = Continuous; 
	}
      // type of jth outcome is set as continuous if any individual has value that is 
      // (not missing) and not (0 or 1)
      // maybe should check for all values missing
      //     if(i == NumIndividuals){
      //       Log << "ERROR: all outcomes missing\n";
      //       exit(1);
      //     }
      
      Log << "Regressing on ";    
      if( OutcomeType[j] == Binary ){
	Log << "Binary variable: ";
	if(numoutcomes==1)RegType = Logistic;
      }
      else if(OutcomeType[j] == Continuous ){
	Log << "Continuous variable: ";
	if(numoutcomes==1)RegType = Linear;
      }
      Log << outcomeVarData_[0][j+Firstcol];
      Log << ".\n";
      OutcomeLabels.push_back(outcomeVarData_[0][j+Firstcol]);
    }
    Log << "\n";
  }
  return RegType;
}

void InputData::CheckCovariatesFile(LogWriter &Log)const{
  if( NumIndividuals != (int)covariatesMatrix_.nRows() - 1 ){
    Log << "ERROR: Genotypes file has " << NumIndividuals << " observations and Covariates file has "
	<< covariatesMatrix_.nRows() - 1 << " observations.\n";
    exit(1);
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

//returns sex value from genotypes file for individual i
Sex InputData::GetSexValue(int i)const{
  //if (options->getgenotypesSexColumn() == 1) {
    int sex = StringConvertor::toInt(geneticData_[i][1]);
    if (sex > 2) {
      cout << "Error: sex must be coded as 0 - missing, 1 - male or 2 - female.\n";
      exit(0);
    }        
    //}
    return (Sex) sex;
}

void InputData::GetGenotype(int i, int SexColumn, const Genome &Loci, vector<genotype>* genotypes, bool** Missing)const{
  unsigned int lociI = 0;
  unsigned complocus = 0;
  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    for(unsigned int j = 0; j < Loci.GetSizeOfChromosome(c); ++j){
      genotype G;
      // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype 
      int numLoci = Loci(complocus)->GetNumberOfLoci();
      
      unsigned int count = 0;
      for (int locus = 0; locus < numLoci; locus++) {
	vector<unsigned short> g(2);
	int col = 1 + SexColumn + lociI;
	if (IsPedFile)col = 1 + SexColumn + 2*lociI;
	  
	StringConvertor::toIntPair(&g, geneticData_[i][col]);

	if(g[0] > Loci(complocus)->GetNumberOfAllelesOfLocus(locus) || (g[1] > Loci(complocus)->GetNumberOfAllelesOfLocus(locus)))
	  throwGenotypeError(i, locus, Loci(complocus)->GetLabel(0), 
			     g[0], g[1], Loci(complocus)->GetNumberOfAllelesOfLocus(locus) );
	lociI++;
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
       << "Individual " << ind << " at locus " << label << locus
       << " has genotype " << g0 << ", " << g1 << " \n"
       << "Number of allelic states at locus = " << numalleles << "\n";
  if(ind == NumIndividuals)
    exit(1);
}

void InputData::getOutcomeTypes(DataType* T)const{
  for(unsigned i = 0; i < outcomeVarMatrix_.nCols(); ++i)
    T[i] = OutcomeType[i];
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

const DataMatrix& InputData::getAlleleFreqMatrix() const
{
    return alleleFreqMatrix_;
}

const DataMatrix& InputData::getHistoricalAlleleFreqMatrix() const
{
    return historicalAlleleFreqMatrix_;
}

const DataMatrix& InputData::getPriorAlleleFreqMatrix() const
{
    return priorAlleleFreqMatrix_;
}

const DataMatrix& InputData::getOutcomeVarMatrix() const
{
    return outcomeVarMatrix_;
}

const DataMatrix& InputData::getReportedAncestryMatrix() const
{
    return reportedAncestryMatrix_;
}

const DataMatrix& InputData::getCovariatesMatrix() const
{
    return covariatesMatrix_;
}
std::string *InputData::GetPopLabels() const{
  return PopulationLabels;
}
Vector_s InputData::getOutcomeLabels()const{
  return OutcomeLabels;
}
std::vector<std::string> InputData::getLocusLabels()const{
  return LocusLabels;
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
  alleleFreqMatrix_.clear();
  historicalAlleleFreqMatrix_.clear();
  priorAlleleFreqMatrix_.clear();
  etaPriorMatrix_.clear();
  reportedAncestryMatrix_.clear();
}
