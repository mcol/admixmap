#include "InputData.h"
#include "AdmixOptions.h"
#include "StringSplitter.h"
#include "StringConvertor.h"

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
        string msg = "Cannt open file for reading: \"";
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

    // If there is no rows, return empty matrix.
    if (0 == numRows) return;

    // Verify that all rows has same length.
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
        ::readFile(options->getGeneInfoFilename(), geneInfoData_);
        ::readFile(options->getGeneticDataFilename(), geneticData_);
        ::readFile(options->getInputFilename(), inputData_);
        ::readFile(options->getTargetFilename(), targetData_);                
        ::readFile(options->getAlleleFreqFilename(), alleleFreqData_);
        ::readFile(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_);            
        ::readFile(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_);
        ::readFile(options->getEtaPriorFilename(), etaPriorData_);
        ::readFile(options->getMLEFilename(), MLEData_);
        ::readFile(options->getReportedAncestryFilename(), reportedAncestryData_);

        // Form matrixes.
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
