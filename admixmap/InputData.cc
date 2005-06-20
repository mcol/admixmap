#include "InputData.h"
#include "AdmixOptions.h"
#include "StringSplitter.h"
#include "StringConvertor.h"
#include "Genome.h"
#include "Chromosome.h"

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

static void readFile(const char *fname, Matrix_s& data)
{
    if (0 == fname || 0 == strlen(fname)) return;

    ifstream in(fname);
    if (!in.is_open()) {
        string msg = "Cannot open file for reading: \"";
        msg += fname;
        msg += "\"";
        throw runtime_error(msg.c_str());
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

void InputData::readData(AdmixOptions *options, LogWriter * /*log*/)
{
  try
    {
      // Read all input files.
      ::readFile(options->getGeneInfoFilename(), geneInfoData_);   //locusfile
      ::readFile(options->getGeneticDataFilename(), geneticData_); //genotypes file
      ::readFile(options->getInputFilename(), inputData_);         //covariates file
      ::readFile(options->getTargetFilename(), targetData_);       //outcomevar file                
      ::readFile(options->getAlleleFreqFilename(), alleleFreqData_);
      ::readFile(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_);            
      ::readFile(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_);
      ::readFile(options->getEtaPriorFilename(), etaPriorData_);
      ::readFile(options->getMLEFilename(), MLEData_);
      ::readFile(options->getReportedAncestryFilename(), reportedAncestryData_);
      
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
  const int isPedFile = 2*geneticData_[0].size() - 1 == geneticData_[1].size() ? 1 : 0;
  options->IsPedFile(isPedFile);

  return (isPedFile==1);
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
void InputData::GetGenotype(int i,AdmixOptions *options,Genome &Loci, unsigned short ****genotype){
  unsigned int lociI = 0;
  
  *genotype = new unsigned short **[Loci.GetNumberOfCompositeLoci()];
  for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); ++j){
    // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype 
    int numLoci = Loci(j)->GetNumberOfLoci();
    
    (*genotype)[j] = new unsigned short *[numLoci];

    for (int locus = 0; locus < numLoci; locus++) {
      (*genotype)[j][locus] = new unsigned short[2];
  
      if (options->IsPedFile() == 1) {
	StringConvertor::toIntPair((*genotype)[j][locus],geneticData_[i][1 + options->getgenotypesSexColumn() + 2*lociI]);
      } 
      else 
	{
	  StringConvertor::toIntPair((*genotype)[j][locus],geneticData_[i][1 + options->getgenotypesSexColumn() + lociI]);
	}
      
      lociI++;
    }
  }
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
