/**
 *   ADMIXMAP
 *   InputData.cc 
 *   Class to read and check all input data files
 *   Copyright (c) 2005 LSHTM
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

/**
 *  Auxilary functions to read data from file.
 */
static bool isWhiteLine(const char *p)
{
    while (*p) {
        if (!isspace(*p++)) {
            return false;
        }
    }

    return true;
}

//Extracts population labels from header line of allelefreq input file
static void getPopLabels(const Vector_s& data, size_t Populations, string **labels)
{
  if(data.size() != Populations+1){cout << "Error in getPopLabels\n";exit(1);}
  *labels = new string[ Populations ];
  StringConvertor S;

    for (size_t i = 0; i < Populations; ++i) {
      (*labels)[i] = S.dequote(data[i+1]);
    }
}
void getLabels(const Vector_s& data, string *labels)
{
  for (size_t i = 0, index = 0; i < data.size(); ++i) {
    labels[index++] = data[i];
  }
}

void InputData::readFile(const char *fname, Matrix_s& data)
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
      Log->logmsg(false,"Loading ");
      Log->logmsg(false,fname);
      Log->logmsg(false,".\n");
    }

    data.clear();
    try {
        StringSplitter splitter;

        string line;        

        while (getline(in, line)) {
            if (!isWhiteLine(line.c_str())) {
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

void InputData::readData(AdmixOptions *options, LogWriter *log)
{
  Log = log;
  try
    {
      // Read all input files.
      readFile(options->getLocusFilename(), locusData_);   //locusfile
      readFile(options->getGenotypesFilename(), geneticData_); //genotypes file
      readFile(options->getCovariatesFilename(), inputData_);         //covariates file
      readFile(options->getOutcomeVarFilename(), outcomeVarData_);       //outcomevar file
      readFile(options->getAlleleFreqFilename(), alleleFreqData_);
      readFile(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_);            
      readFile(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_);
      readFile(options->getEtaPriorFilename(), etaPriorData_);
      readFile(options->getMLEFilename(), MLEData_);
      readFile(options->getReportedAncestryFilename(), reportedAncestryData_);

      Log->logmsg(false,"\n");      
      // Form matrices.
      convertMatrix(locusData_, locusMatrix_);
      ::convertMatrix(outcomeVarData_, outcomeVarMatrix_);
      ::convertMatrix(inputData_,  covariatesMatrix_);
      ::convertMatrix(alleleFreqData_, alleleFreqMatrix_);
      ::convertMatrix(historicalAlleleFreqData_, historicalAlleleFreqMatrix_);
      ::convertMatrix(priorAlleleFreqData_, priorAlleleFreqMatrix_);
      ::convertMatrix(etaPriorData_, etaPriorMatrix_);
      ::convertMatrix(MLEData_, MLEMatrix_);
      ::convertMatrix(reportedAncestryData_, reportedAncestryMatrix_);
      
    } catch (const exception& e) {
    cerr << "Exception occured during parsing input file: \n" << e.what() << endl;
    exit(1);
  }
  NumSimpleLoci = getNumberOfSimpleLoci();
  NumCompositeLoci = determineNumberOfCompositeLoci();
  NumIndividuals = getNumberOfIndividuals();

  IsPedFile = determineIfPedFile( options );
  CheckGeneticData(options);
  checkLocusFile(options->getgenotypesSexColumn());
  locusMatrix_ = locusMatrix_.SubMatrix(1, locusMatrix_.nRows() - 1, 1, 2);//remove header and first column of locus file
  if ( strlen( options->getOutcomeVarFilename() ) != 0 )
    options->setRegType( CheckOutcomeVarFile( options->getNumberOfOutcomes(), options->getTargetIndicator()) );
  if ( strlen( options->getCovariatesFilename() ) != 0 )
    CheckCovariatesFile();//detects regression model
  if ( strlen( options->getReportedAncestryFilename() ) != 0 )
    CheckRepAncestryFile(options->getPopulations());
  CheckAlleleFreqs(options);
  
  if(NumIndividuals > 1){
    Log->logmsg(false, NumIndividuals);Log->logmsg(false, " individuals\n");
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
     if( locusMatrix_.get( i, 2 ) == 0.0 ) NumberOfCompositeLoci--;
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
      cerr << "Wrong number of entries in line "<<i+1<<" of genotypesfile" << endl;
      exit(1);
    }
    
  }
}

void InputData::checkLocusFile(int sexColumn)const{
  // Check that loci labels in locusfile are unique and that they match the names in the genotypes file.
  
  for (size_t i = 1; i < locusData_.size(); ++i) {//rows of locusfile
    //check distances are not negative
    if(locusMatrix_.get(i,2) < 0.0){
      cerr<<"Error: distance on line "<<i<<" of locusfile is negative."<<endl;
      exit(1);
    }
    // Check loci names are unique    
    for (size_t j = i + 1; j < locusData_.size(); ++j) {   
      if (locusData_[i][0] == locusData_[j][0]) {
	cerr << "Error in locusfile. Two different loci have the same name. "
	     << locusData_[i][0] << endl;
	exit(2);            
      }
    }
  }

  const size_t numLoci = locusData_.size() - 1;//number of simple loci

  // Compare loci names in locus file and genotypes file.
  for (size_t i = 1; i <= numLoci; ++i) {
    if (locusData_[i][0] != geneticData_[0][i + sexColumn]) {
      cout << "Error. Loci names in locus file and genotypes file are not the same." << endl;
      cout << "Loci names causing an error are: " << locusData_[i][0] << " and " 
	   << geneticData_[0][i + sexColumn] << endl;
      //cout << options->getgenotypesSexColumn() << endl;
      exit(2);
    }
  } 
}

//checks consistency of supplied allelefreqs with locusfile
//and determines number of populations and population labels
void InputData::CheckAlleleFreqs(AdmixOptions *options){
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
    while( index < locusMatrix_.nRows() - 1 && locusMatrix_.get( index, 1 ) == 0 );
    NumberOfStates += states;
  }


  //fixed allele freqs
  if( strlen( options->getAlleleFreqFilename() ) ){
    freqtype = "";
    infile = true;
    nrows = alleleFreqMatrix_.nRows()-1;
    expectednrows = NumberOfStates-NumCompositeLoci;
    Populations = alleleFreqMatrix_.nCols() - 1;// -1 for ids in first col
    ::getPopLabels(alleleFreqData_[0], Populations, &PopulationLabels);
  }
  
  //Historic allelefreqs
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    freqtype = "historic";
    infile = true;
    nrows = historicalAlleleFreqMatrix_.nRows();
    expectednrows = NumberOfStates+1;
    Populations = historicalAlleleFreqMatrix_.nCols() - 1;
    ::getPopLabels(historicalAlleleFreqData_[0], Populations, &PopulationLabels);

  }
  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
      freqtype = "prior";
      infile = true;
      nrows = priorAlleleFreqMatrix_.nRows();
      expectednrows = NumberOfStates+1;
      Populations = priorAlleleFreqMatrix_.nCols() - 1;
      ::getPopLabels(priorAlleleFreqData_[0], Populations, &PopulationLabels);
  }
  if(infile){
    if(nrows != expectednrows){
      Log->logmsg(true,"Incorrect number of rows in ");
      Log->logmsg(true, freqtype);
      Log->logmsg(true, "allelefreqfile.\n");
      Log->logmsg(true,"Expecting ");
      Log->logmsg(true,expectednrows);
      Log->logmsg(true," rows, but there are ");
      Log->logmsg(true,nrows);
      Log->logmsg(true," rows.\n");
      exit(1);
    }
    options->setPopulations(Populations);
  }
  else{//'populations' option
    if(Populations < 1){
      Log->logmsg(true, "ERROR: populations = ");
      Log->logmsg(true, options->getPopulations());
      Log->logmsg(true, "\n");
      exit(1);
    }
    PopulationLabels = new string[ Populations ];
    for( int j = 0; j < Populations; j++ ){
      stringstream poplabel;
      string result;
      poplabel << "Pop" << j+1;
      result = poplabel.str();
      PopulationLabels[j] = result;
    }
//     for( int i = 0; i < NumberOfCompositeLoci; i++ ){
//       if(Loci->GetNumberOfStates(i) < 2){
// 	Log->logmsg(true, "ERROR: The number of alleles at a locus is ");
// 	Log->logmsg(true, Loci->GetNumberOfStates(i));
// 	Log->logmsg(true, "\n");
// 	exit(1);
//       }
//     }
  }
}

RegressionType InputData::CheckOutcomeVarFile(int NumOutcomes, int Firstcol){
  //check outcomevarfile and genotypes file have the same number of cols
  if( (int)outcomeVarMatrix_.nRows() - 1 != NumIndividuals ){
    Log->logmsg(true,"ERROR: Genotypes file has ");
    Log->logmsg(true,NumIndividuals);
    Log->logmsg(true," observations and Outcomevar file has ");
    Log->logmsg(true,outcomeVarMatrix_.nRows() - 1);
    Log->logmsg(true," observations.\n");
    exit(1);
  }
  //check the number of outcomes specified is not more than the number of cols in outcomevarfile
  int numoutcomes = NumOutcomes;
  if(NumOutcomes > -1){//options 'numberofregressions' used
    if((int)outcomeVarMatrix_.nCols() - Firstcol < NumOutcomes){
      numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;//adjusts if too large
      Log->logmsg(true, "ERROR: 'outcomes' is too large, setting to ");
      Log->logmsg(true, numoutcomes);
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
      
      for(int i = 0; i < NumIndividuals; ++i)
	if(!outcomeVarMatrix_.isMissing(i, j) && (outcomeVarMatrix_.get( i, j ) == 0 || outcomeVarMatrix_.get( i, j ) == 1) )
	  OutcomeType[j] = Binary;
	else OutcomeType[j] = Continuous;
      //in this way, the outcome type is set as binary only if all individuals have outcome values of 1 or 0
      //otherwise, a continuous outcome of 1.0 or 0.0 could lead to the type being wrongly set to binary.
      
      //need to check for allmissing
      //     if(i == NumIndividuals){
      //       Log->logmsg(true, "ERROR: all outcomes missing\n");
      //       exit(1);
      //     }
      
      Log->logmsg(true,"Regressing on ");    
      if( OutcomeType[j] == Binary ){
	Log->logmsg(true,"Binary variable: ");
	if(numoutcomes==1)RegType = Logistic;
      }
      else if(OutcomeType[j] == Continuous ){
	Log->logmsg(true,"Continuous variable: ");
	if(numoutcomes==1)RegType = Linear;
      }
      Log->logmsg(true,outcomeVarData_[0][j+Firstcol]);
      Log->logmsg(true,".\n");
      OutcomeLabels.push_back(outcomeVarData_[0][j+Firstcol]);
    }
    Log->logmsg(true, "\n");
  }
  return RegType;
}

void InputData::CheckCovariatesFile()const{
  if( NumIndividuals != (int)covariatesMatrix_.nRows() - 1 ){
    Log->logmsg(true,"ERROR: Genotypes file has ");
    Log->logmsg(true,NumIndividuals);
    Log->logmsg(true," observations and Covariates file has ");
    Log->logmsg(true,covariatesMatrix_.nRows() - 1);
    Log->logmsg(true," observations.\n");
    exit(1);
  }
}

void InputData::CheckRepAncestryFile(int populations)const{
  if( (int)reportedAncestryMatrix_.nRows() != 2 * NumIndividuals ){
    Log->logmsg(false,"ERROR: ");
    Log->logmsg(false,"ReportedAncestry file");
    Log->logmsg(false," has ");
    Log->logmsg(false,reportedAncestryMatrix_.nRows());
    Log->logmsg(false," rows\n");
    Log->logmsg(false,"Genotypesfile");
    Log->logmsg(false," has ");
    Log->logmsg(false,NumIndividuals);
    Log->logmsg(false," rows\n");
    exit(1);}
  if( (int)reportedAncestryMatrix_.nCols() != populations ){
    Log->logmsg(false,"ERROR: ");
    Log->logmsg(false,"ReportedAncestry file");
    Log->logmsg(false," has ");
    Log->logmsg(false,reportedAncestryMatrix_.nCols());
    Log->logmsg(false," cols\n");
    Log->logmsg(false, "AlleleFreq file");
    Log->logmsg(false," has ");
    Log->logmsg(false,populations);
    Log->logmsg(false," cols\n");
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

void InputData::GetGenotype(int i, int SexColumn, const Genome &Loci, unsigned short ****genotype)const{
  unsigned int lociI = 0;
  
  *genotype = new unsigned short **[Loci.GetNumberOfCompositeLoci()];

    for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); ++j){
      // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype 
      int numLoci = Loci(j)->GetNumberOfLoci();
      
      (*genotype)[j] = new unsigned short *[numLoci];
      
      for (int locus = 0; locus < numLoci; locus++) {
	(*genotype)[j][locus] = new unsigned short[2];
	int col = 1 + SexColumn + lociI;
	if (IsPedFile)col = 1 + SexColumn + 2*lociI;
	  
	StringConvertor::toIntPair((*genotype)[j][locus],geneticData_[i][col]);

	if((*genotype)[j][locus][0] > Loci(j)->GetNumberOfAllelesOfLocus(locus) ||
	   (*genotype)[j][locus][1] > Loci(j)->GetNumberOfAllelesOfLocus(locus))
	  throwGenotypeError(i, locus, Loci(j)->GetLabel(j), 
			     (*genotype)[j][locus][0], (*genotype)[j][locus][1], Loci(j)->GetNumberOfAllelesOfLocus(locus) );
	lociI++;
      }
    }
}

void InputData::throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles)const{
  Log->logmsg(false, "Error in genotypes file:\n");
  Log->logmsg(false, "Individual ");
  Log->logmsg(false, ind);
  Log->logmsg(false, " at locus ");
  Log->logmsg(false, label);Log->logmsg(false, locus);
  Log->logmsg(false, " has genotype ");
  Log->logmsg(false, g0);Log->logmsg(false, ", ");
  Log->logmsg(false, g1);Log->logmsg(false, " \n");
  Log->logmsg(false, "Number of allelic states at locus = ");
  Log->logmsg(false, numalleles);Log->logmsg(false, "\n");
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

const Matrix_s& InputData::getMLEData() const
{
    return MLEData_;
}

const Matrix_s& InputData::getReportedAncestryData() const
{
    return reportedAncestryData_;
}

const DataMatrix& InputData::getEtaPriorMatrix() const
{
    return etaPriorMatrix_;
}

const DataMatrix& InputData::getMLEMatrix() const
{
    return MLEMatrix_;
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
