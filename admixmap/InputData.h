// *-*-C++-*-*
#ifndef INPUT_DATA_H
#define INPUT_DATA_H 1

#include "common.h"

/**
 *  Forward declarations.
 */    
class AdmixOptions;
class LogWriter;


class InputData 
{
public:    

    /**
     *  Constructor.
     */    
    InputData();
    
    /**
     *  Descructor.
     */    
    ~InputData();
    

    /**
     *  Read all input data and store in internal structures.
     */    
    void readData(AdmixOptions *options, LogWriter *log);    

    /**
     *  Getters to retrieve data (in string form).
     */
    const Matrix_s& getGeneInfoData() const;
    const Matrix_s& getGeneticData()  const;
    const Matrix_s& getInputData()    const;
    const Matrix_s& getTargetData()   const;
    const Matrix_s& getAlleleFreqData() const;
    const Matrix_s& getHistoricalAlleleFreqData() const;
    const Matrix_s& getPriorAlleleFreqData() const;
    const Matrix_s& getEtaPriorData() const;
    const Matrix_s& getMLEData() const;
    const Matrix_s& getReportedAncestryData() const;

    /**
     *  Getters to retrieve data (converted to Matrix_d).
     */    
    const Matrix_d& getGeneInfoMatrix() const;
    const Matrix_d& getTargetMatrix() const;
    const Matrix_d& getInputMatrix() const;
    const Matrix_d& getAlleleFreqMatrix() const;
    const Matrix_d& getHistoricalAlleleFreqMatrix() const;
    const Matrix_d& getPriorAlleleFreqMatrix() const;
    const Matrix_d& getEtaPriorMatrix() const;
    const Matrix_d& getMLEMatrix() const;
    const Matrix_d& getReportedAncestryMatrix() const;

private:    
    Matrix_s geneInfoData_;
    Matrix_s geneticData_;
    Matrix_s inputData_;
    Matrix_s targetData_;
    Matrix_s alleleFreqData_;
    Matrix_s historicalAlleleFreqData_;
    Matrix_s priorAlleleFreqData_;
    Matrix_s etaPriorData_;
    Matrix_s MLEData_;
    Matrix_s reportedAncestryData_;

    Matrix_d geneInfoMatrix_;
    Matrix_d inputMatrix_;
    Matrix_d targetMatrix_;
    Matrix_d alleleFreqMatrix_;
    Matrix_d historicalAlleleFreqMatrix_;
    Matrix_d priorAlleleFreqMatrix_;
    Matrix_d etaPriorMatrix_;
    Matrix_d MLEMatrix_;
    Matrix_d reportedAncestryMatrix_;

private:
    /**
     *  UNIMPLEMENTED: to avoid undesired copying.
     */    
    InputData(const InputData&);
    void operator=(const InputData&);
};

#endif /* !defined INPUT_DATA_H */
