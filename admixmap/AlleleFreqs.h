// *-*-C++-*-*
#ifndef ALLELEFREQS_H
#define ALLELEFREQS_H 1
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"
#include "LocusVisitor.h"
#include "AlleleFreqOutputter.h"

class AlleleFreqs{

public:
  AlleleFreqs();
  ~AlleleFreqs();
  void Initialise(AdmixOptions *options, const Matrix_d& etaprior,LogWriter *Log,std::string *PopulationLabels, double rho);
  void load_f(double rho,Genome *chrm);
  void UpdateAlleleFreqs(int iteration,int);
  void InitializeOutputFile(AdmixOptions *options, std::string *PopulationLabels);
//outputs SumEta  to ErgodicAverageFile
  void OutputErgodicAvg( int iteration,AdmixOptions *options, std::ofstream *avgstream);
  //outputs eta to dispparamfile
  void OutputEta(int iteration, AdmixOptions *options, std::ofstream *LogFileStreamPtr);

  Genome *getLoci();

  int GetNumberOfCompositeLoci();

  void accept();

  void OutputFST(bool IsPedFile);

  void LoadAlleleFreqs(AdmixOptions *options, Genome **chrm,LogWriter *Log, InputData *data,std::string **PopulationLabels);

  void Reset();//resets Loci object

  Vector_d *geteta();
  Vector_d *getSumEta();
  Vector_d getLociCorrSummary();

private:
  int Number, Populations;
  Vector_d eta; //dispersion parameter
  Vector_d psi,tau;//parameters of gamma prior on eta 
  /**
   * eta has Gamma prior with shape and scale parameters phi and tau
   */

  Genome Loci;
  TuneRW *TuneEtaSampler;
  int w;
  Vector_i NumberAccepted;
  Vector_d etastep;
  double etastep0;

  Vector_d SumEta;
  Vector_d SumAcceptanceProb;
  double psi0;

  // LociCorrSummary is a vector of terms of the form exp( - rho*x_i) where x_i is the map distance between two adjacent loci
  // with a global rho model, this vector is same for all individuals and calculated only once. 
  Vector_d LociCorrSummary;//summary of correlation in ancestry between loci,was called f in Latent

  bool isHistoricAlleleFreq;//indicator for dispersion model
  LocusVisitor* allelefreqoutput;// object to output allele frequencies
  std::ofstream outputstream;//outputs eta to paramfile
  std::ofstream fstoutputstream;

  void checkLociNames(AdmixOptions *,InputData *data_);
    
  void getLabels( const std::string, Vector_i, std::string* );
 
  void loadAlleleStatesAndDistances(std::vector<std::string>* ChrmLabels,AdmixOptions *,InputData *data_,LogWriter *);

  void OpenFSTFile(AdmixOptions *options,LogWriter *Log); 

  static double strangExp( double );  
};






#endif
