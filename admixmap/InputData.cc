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
 *  Auxilary function that converts Matrix_s to Matrix_d
 */
static void convertMatrix(const Matrix_s& data, Matrix_d& m)
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
    m.SetNumberOfElements(numRows, numCols);
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            if (StringConvertor::isMissingValue(data[i][j])) {
                m.SetMissingElement(i, j);
            } else {
                m(i, j) = StringConvertor::toFloat(data[i][j]);
            }
        }
    }
}


/**
 *  InputData members.
 */

InputData::InputData()
{
}

InputData::~InputData()
{
}

void InputData::readData(AdmixOptions *options, LogWriter *log)
{
  Log = log;
  try
    {
      // Read all input files.
      readFile(options->getGeneInfoFilename(), geneInfoData_);   //locusfile
      readFile(options->getGeneticDataFilename(), geneticData_); //genotypes file
      readFile(options->getInputFilename(), inputData_);         //covariates file
      readFile(options->getTargetFilename(), targetData_);       //outcomevar file                
      readFile(options->getAlleleFreqFilename(), alleleFreqData_);
      readFile(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_);            
      readFile(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_);
      readFile(options->getEtaPriorFilename(), etaPriorData_);
      readFile(options->getMLEFilename(), MLEData_);
      readFile(options->getReportedAncestryFilename(), reportedAncestryData_);
      
      // Form matrices.
      convertMatrix(geneInfoData_, geneInfoMatrix_);
      if (options->getTextIndicator()) {
	geneInfoMatrix_.SubMatrix2(1, geneInfoMatrix_.GetNumberOfRows() - 1, 1, 2);
      } else {
	geneInfoMatrix_.SubMatrix2(0, geneInfoMatrix_.GetNumberOfRows() - 1, 0, 1);
      }
      
      ::convertMatrix(targetData_, targetMatrix_);
      ::convertMatrix(inputData_,  inputMatrix_);
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
  NumIndividuals = getNumberOfIndividuals();
  IsPedFile = determineIfPedFile( options );
  CheckGeneticData(options->getgenotypesSexColumn());
  checkLociNames(options);
 if ( strlen( options->getTargetFilename() ) != 0 )
   CheckOutcomeVarFile((bool)(options->getAnalysisTypeIndicator() != 5));
 if ( strlen( options->getInputFilename() ) != 0 )
   CheckCovariatesFile();
 if ( strlen( options->getReportedAncestryFilename() ) != 0 )
   CheckRepAncestryFile(options->getPopulations());
  //convertGenotypesToIntArray(options );
}

//determine number of individuals by counting lines in genotypesfile 
int InputData::getNumberOfIndividuals() {
  return(geneticData_.size() - 1);
}

//determine number of loci by counting rows of locusfile
int InputData::getNumberOfSimpleLoci() {
  return(geneInfoData_.size() - 1);
}
bool InputData::determineIfPedFile(AdmixOptions *options) {
  // Determine if genotype table is in pedfile format by testing if number of strings in row 1 equals
  // twice the number of strings in the header row minus one. 
  // 
  const bool isPedFile = (bool)(2*geneticData_[0].size() - 1 == geneticData_[1].size());
  options->IsPedFile(isPedFile);

  return (isPedFile);
}

//checks number of loci in genotypes file is the same as in locusfile
void InputData::CheckGeneticData(int genotypesSexColumn){
  for(int i = 1; i <= NumIndividuals; ++i){
    //should use logmsg
    if (IsPedFile) {
      if ((int)geneticData_[i].size()-1 != 2*NumSimpleLoci + genotypesSexColumn) {
	cout << "Error in formatting of line " <<i+1<<" of genotypesfile"<< endl;
	exit(0);
      }
    } else {
      if ((int)geneticData_[i].size()-1 != NumSimpleLoci + genotypesSexColumn) {
	cout << "Error in formatting of line "<<i+1<<" of genotypesfile" << endl;
	exit(0);
      }
    }
  }
}

void InputData::checkLociNames(AdmixOptions *options){
  // Check that loci labels in locusfile are unique and that they match the names in the genotypes file.
  
  // Check loci names are unique    
  for (size_t i = 1; i < geneInfoData_.size(); ++i) {
    for (size_t j = i + 1; j < geneInfoData_.size(); ++j) {   
      if (geneInfoData_[i][0] == geneInfoData_[j][0]) {
	cerr << "Error in locusfile. Two different loci have the same name. "
	     << geneInfoData_[i][0] << endl;
	exit(2);            
      }
    }
  }

  const size_t numLoci = geneInfoData_.size() - 1;

//this should be in InputData
    // Determine if "Sex" column present in genotypes file.
    if (numLoci == geneticData_[0].size() - 1) {
        options->setgenotypesSexColumn(0);
    } else if (numLoci == geneticData_[0].size() - 2) {
        options->setgenotypesSexColumn(1);
    } else {
        cerr << "Error. Number of loci in genotypes file does not match number in locus file." << endl;
        exit(2);
    }

    // Compare loci names in locus file and genotypes file.
    for (size_t i = 1; i <= numLoci; ++i) {
        if (geneInfoData_[i][0] != geneticData_[0][i + options->getgenotypesSexColumn()]) {
            cout << "Error. Loci names in locus file and genotypes file are not the same." << endl;
            cout << "Loci names causing an error are: " << geneInfoData_[i][0] << " and " 
                 << geneticData_[0][i + options->getgenotypesSexColumn()] << endl;
            cout << options->getgenotypesSexColumn() << endl;
            exit(2);
        }
    } 
}

//checks consistency of supplied allelefreqs with locusfile
void InputData::CheckAlleleFreqs(AdmixOptions *options, int NumberOfCompositeLoci, int NumberOfStates){
  string freqtype = "x";
  bool infile = false;
  int nrows, expectednrows;

  //fixed allele freqs
  if( strlen( options->getAlleleFreqFilename() ) ){
    freqtype = "";
    infile = true;
    nrows = alleleFreqMatrix_.GetNumberOfRows()-1;
    expectednrows = NumberOfStates-NumberOfCompositeLoci;
  }
  
  //Historic AlleleFreqs
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    freqtype = "historic";
    infile = true;
    nrows = historicalAlleleFreqMatrix_.GetNumberOfRows();
    expectednrows = NumberOfStates+1;
  }
  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
      freqtype = "prior";
      infile = true;
      nrows = priorAlleleFreqMatrix_.GetNumberOfRows();
      expectednrows = NumberOfStates+1;
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
  }
  else{//'populations' option
    if(options->getPopulations() < 1){
      Log->logmsg(true, "ERROR: populations = ");
      Log->logmsg(true, options->getPopulations());
      Log->logmsg(true, "\n");
      exit(1);
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

void InputData::CheckOutcomeVarFile(bool singleRegression){
  if( targetMatrix_.GetNumberOfRows() - 1 != NumIndividuals ){
    Log->logmsg(true,"ERROR: Genotypes file has ");
    Log->logmsg(true,NumIndividuals);
    Log->logmsg(true," observations and Outcomevar file has ");
    Log->logmsg(true,inputMatrix_.GetNumberOfRows() - 1);
    Log->logmsg(true," observations.\n");
    exit(1);
  }
  if(singleRegression){
    if( targetMatrix_.GetNumberOfRows() - 1 != NumIndividuals ){
      Log->logmsg(true,"Outcomevar file has ");
      Log->logmsg(true,targetMatrix_.GetNumberOfRows() - 1);
      Log->logmsg(true," observations and Genotypes file has ");
      Log->logmsg(true,NumIndividuals);
      Log->logmsg(true," observations.\n");
      exit(1);
    }
  }
}

void InputData::CheckCovariatesFile(){
  if( NumIndividuals != inputMatrix_.GetNumberOfRows() - 1 ){
    Log->logmsg(true,"ERROR: Genotypes file has ");
    Log->logmsg(true,NumIndividuals);
    Log->logmsg(true," observations and Covariates file has ");
    Log->logmsg(true,inputMatrix_.GetNumberOfRows() - 1);
    Log->logmsg(true," observations.\n");
    exit(1);
  }
}

void InputData::CheckRepAncestryFile(int populations){
  if( reportedAncestryMatrix_.GetNumberOfRows() != 2 * NumIndividuals ){
    Log->logmsg(false,"ERROR: ");
    Log->logmsg(false,"ReportedAncestry file");
    Log->logmsg(false," has ");
    Log->logmsg(false,reportedAncestryMatrix_.GetNumberOfRows());
    Log->logmsg(false," rows\n");
    Log->logmsg(false,"Genotypesfile");
    Log->logmsg(false," has ");
    Log->logmsg(false,NumIndividuals);
    Log->logmsg(false," rows\n");
    exit(1);}
  if( reportedAncestryMatrix_.GetNumberOfCols() != populations ){
    Log->logmsg(false,"ERROR: ");
    Log->logmsg(false,"ReportedAncestry file");
    Log->logmsg(false," has ");
    Log->logmsg(false,reportedAncestryMatrix_.GetNumberOfCols());
    Log->logmsg(false," cols\n");
    Log->logmsg(false, "AlleleFreq file");
    Log->logmsg(false," has ");
    Log->logmsg(false,populations);
    Log->logmsg(false," cols\n");
    exit(1);
  }
}

//returns sex value from genotypes file for individual i
int InputData::GetSexValue(int i){
  //if (options->getgenotypesSexColumn() == 1) {
    int sex = StringConvertor::toInt(geneticData_[i][1]);
    if (sex > 2) {
      cout << "Error: sex must be coded as 0 - missing, 1 - male or 2 - female.\n";
      exit(0);
    }        
    //}
    return sex;
}


//converts genotypes stored as Matrix_s strings to genotypes stored as Matrix g integer pairs
//also removes cols for ID and sex - this should have been a separate step  
void InputData::convertGenotypesToIntArray(AdmixOptions *options ) {
  int firstcol = 1 + options->getgenotypesSexColumn();
  genotype g;
  Vector_g vgenotypes; // vector of individual's genotypes
  unsigned int *a = new unsigned int[2];

  //loop over individuals
  for( unsigned int indiv=0; indiv < geneticData_.size() - 1; indiv++ ) {
    Vector_s & genotypes_s = geneticData_[indiv + 1];
    // loop over simple loci to store genotype strings as unsigned integers 
    for( int locus = 0; locus < NumSimpleLoci; locus++ ) {
      // needs fixing to deal with X chr data in males
      // should throw exception for index out of range 
      StringConvertor::toIntPair(a, genotypes_s[firstcol + locus].c_str() );
      g.alleles[0] = a[0];
      g.alleles[1] = a[1];
    }
    // append genotype to vgenotypes
    vgenotypes.push_back( g );
  }
  // append genotypes vector to genotypes_g
  genotypes_g.push_back( vgenotypes );
  delete[] a;
}


// loop over simple loci within each composite locus to store genotypes at S simple loci within 
// each composite locus as vector of length 2S.  
void InputData::convertToVectorsOverCLoci(Genome & Loci, Chromosome **chrm) {
  int NumChromosomes = Loci.GetNumberOfChromosomes();
  int NumCLoci = Loci.GetNumberOfCompositeLoci();
  int NumSimpleLoci;
  int CLocus;
  std::vector< std::vector<unsigned int> > genotypes_cloci; 
  for( int indiv=0; indiv < NumIndividuals; indiv++ ) { // loop over individuals
    for( int j = 0; j < NumChromosomes; j++ ){ // loop over chromosomes
      // loop over composite loci
      for( int jj = 0; jj < NumCLoci; jj++ ){
     	CLocus = chrm[j]->GetLocus(jj); // get number of this composite locus
     	NumSimpleLoci = Loci(CLocus)->GetNumberOfLoci(); // 
     	// create new vector genotypes.clocus for this composite locus 
     	std::vector<unsigned int> genotypes_clocus(NumSimpleLoci * 2, 0); // should delete after appending to genotypes_cloci
     	// loop over simple loci within composite locus to assign elements of vector genotypes.clocus
     	for (int locus=0; locus < NumSimpleLoci; locus++) { //
     	  genotypes_clocus[locus*2]    = genotypes_g[indiv][locus].alleles[0];
     	  genotypes_clocus[locus*2+1]  = genotypes_g[indiv][locus].alleles[1];
     	}
     	// append genotypes for composite locus to vector over composite loci for this individual
	genotypes_cloci.push_back( genotypes_clocus);
	genotypes_clocus.clear(); // do not re-use this object as it would have to be resized - possible memory leaks
      }
      // append genotypes for individual to vector over individuals
      genotypes_c.push_back(genotypes_cloci);
      genotypes_cloci.clear(); // could re-use this object
    }
  }
}

//TODO: maybe have numChromosomes, NumLoci etc members of InputData
void InputData::GetGenotype(int i, int SexColumn, Genome &Loci, unsigned short ****genotype){
  unsigned int lociI = 0;
  
  *genotype = new unsigned short **[Loci.GetNumberOfCompositeLoci()];

    for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); ++j){
      // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype 
      int numLoci = Loci(j)->GetNumberOfLoci();
      
      (*genotype)[j] = new unsigned short *[numLoci];
      
      for (int locus = 0; locus < numLoci; locus++) {
	(*genotype)[j][locus] = new unsigned short[2];
	
	if (IsPedFile) {
	  StringConvertor::toIntPair((*genotype)[j][locus],geneticData_[i][1 + SexColumn + 2*lociI]);
	} 
	else 
	  {
	    StringConvertor::toIntPair((*genotype)[j][locus],geneticData_[i][1 + SexColumn + lociI]);
	  }
	if((*genotype)[j][locus][0] > Loci(j)->GetNumberOfAllelesOfLocus(locus) ||
	   (*genotype)[j][locus][1] > Loci(j)->GetNumberOfAllelesOfLocus(locus))
	  throwGenotypeError(i, locus, Loci(j)->GetLabel(j), 
			     (*genotype)[j][locus][0], (*genotype)[j][locus][1], Loci(j)->GetNumberOfAllelesOfLocus(locus) );
	lociI++;
      }
    }
}

void InputData::throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles){
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

const Matrix_s& InputData::getGeneInfoData() const
{
    return geneInfoData_;
}

const Matrix_s& InputData::getGeneticData() const
{
    return geneticData_;
}

const Matrix_s& InputData::getInputData() const
{
    return inputData_;
}

const Matrix_s& InputData::getTargetData() const
{
    return targetData_;
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

const Matrix_d& InputData::getEtaPriorMatrix() const
{
    return etaPriorMatrix_;
}

const Matrix_d& InputData::getMLEMatrix() const
{
    return MLEMatrix_;
}

const Matrix_d& InputData::getGeneInfoMatrix() const
{
    return geneInfoMatrix_;
}

const Matrix_d& InputData::getAlleleFreqMatrix() const
{
    return alleleFreqMatrix_;
}

const Matrix_d& InputData::getHistoricalAlleleFreqMatrix() const
{
    return historicalAlleleFreqMatrix_;
}

const Matrix_d& InputData::getPriorAlleleFreqMatrix() const
{
    return priorAlleleFreqMatrix_;
}

const Matrix_d& InputData::getTargetMatrix() const
{
    return targetMatrix_;
}

const Matrix_d& InputData::getReportedAncestryMatrix() const
{
    return reportedAncestryMatrix_;
}

const Matrix_d& InputData::getInputMatrix() const
{
    return inputMatrix_;
}
