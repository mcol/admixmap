// *-*-C++-*-*
#ifndef ALLELEFREQS_H
#define ALLELEFREQS_H 1
#include "InputData.h"
//#include "Chromosome.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"
#include "LocusVisitor.h"
#include "AlleleFreqOutputter.h"


class AlleleFreqs{
  
public:
  AlleleFreqs();
  ~AlleleFreqs();
  void Initialise(AdmixOptions *options, const Matrix_d& etaprior, LogWriter *Log, 
		  std::string *PopulationLabels, double rho);
  
  void load_f(double rho,Chromosome **chrm); // should be moved to chromosome object
  
  void Update(int iteration,int);

  // initialize output file for samples of  
  // void InitializeOutputFile(AdmixOptions *options, std::string *PopulationLabels);
  
  //initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log);

  //outputs ergodic averages of dispersion parameters
  void OutputErgodicAvg( int iteration,AdmixOptions *options, std::ofstream *avgstream);
  
  //output samples of dispersion parameters to dispparamfile
  void OutputEta(int iteration, AdmixOptions *options, std::ofstream *LogFileStreamPtr);
  
  Genome *getLoci();
  CompositeLocus *getLocus(int);
  
  int GetNumberOfCompositeLoci();
  
  void OutputAlleleFreqs();
  
  void OutputFST(bool IsPedFile);
  
  void LoadAlleleFreqs(AdmixOptions *options, Chromosome ***chrm,LogWriter *Log, InputData *data,std::string **PopulationLabels);
  
  void Reset();  // resets Loci object
  int IsRandom(); // why not bool?
  Vector_d GetFst(int locus);
  void UpdateFst();
  Vector_d GetStatsForEta( int , int locus);
  double GetAlleleProbsMAP( int x, int ancestry , int locus);
  Vector_d GetPriorAlleleFreqs( int locus, int population );
  

  Vector_i GetAlleleCounts( int locus, int population );

  Vector_d getAlleleFreqsMAP( int locus, int population );

  Matrix_d &GetAlleleFreqs(int locus);

  Matrix_i &GetAlleleCounts(int locus);

  
  // does not need a public method - this belongs in CompositeLocus object
  int GetNumberOfStates(int locus);

  void UpdateAlleleCounts(int locus, Vector_i Haplotypes, Vector_i ancestry );

  void UpdateAlleleCounts_HaploidData(int locus, const vector<unsigned int>& genotype, int ancestry );
  void ResetSumAlleleFreqs();
  void setAlleleFreqsMAP();
  Matrix_d GetLikelihood( int locus, const vector<unsigned int> genotype, Vector_i Haplotypes, bool diploid, bool fixed);
  Vector_d *geteta();
  Vector_d *getSumEta();

  Vector_d getLociCorrSummary(); // should be in Genome object

  // function to merge rare haplotypes for construction of score tests
  void SetMergedHaplotypes(Vector_d *alpha0, std::ofstream *LogFileStreamPtr, bool IsPedFile);
  
private:
  int Number, Populations;
  Vector_d eta;  //dispersion parameter
  Vector_d psi,tau; // eta has Gamma prior with shape and scale parameters psi and tau
  
  MatrixArray_d Freqs; // allele frequencies except for last allele
  MatrixArray_d AlleleFreqsMAP; // posterior mode of allele freqs
  MatrixArray_d HistoricAlleleFreqs; 
  MatrixArray_d AlleleProbs; // allele freqs including last allele
  
  MatrixArray_i AlleleCounts;
  MatrixArray_d HistoricLikelihoodAlleleFreqs;
  MatrixArray_d PriorAlleleFreqs;

  MatrixArray_d SumAlleleFreqs; // used to compute ergodic average

  // next function can be private as no longer called by AlleleFreqOutputter
  Matrix_d AlleleFreqs::GetSumAlleleFreqs(int locus);
  //Matrix_d SumEta;
  
  Matrix_d Fst;
  Matrix_d SumFst;
  bool IsHistoricAlleleFreq;
  int RandomAlleleFreqs; //indicator for whether allele freqs are fixed or random - should be bool?
 
  Genome Loci; // is this where the Loci object is instantiated?
  
  TuneRW *TuneEtaSampler;
  int w; // the eta sampler is tuned every w updates
  Vector_i NumberAccepted;
  Vector_d etastep;
  double etastep0;
  
  Vector_d SumEta;
  Vector_d SumAcceptanceProb;
  double psi0;
  Vector_d pp; //used to set merged haplotypes, which are used in the allelic association test
  
  //    DARS SampleMu;
  std::vector<TuneRW> *MuProposal;
  
  // LociCorrSummary is a vector of terms of the form exp( - rho*x_i) where x_i is 
  // map distance between two adjacent loci
  // with a global rho model, this vector is same for all individuals and calculated only once. 
  // should be in Genome object
  Vector_d LociCorrSummary; //summary of correlation in ancestry between loci,was called f in Latent
  
  bool isHistoricAlleleFreq; //indicator for dispersion model
  LocusVisitor* allelefreqoutput; // object to output allele frequencies
  std::ofstream outputstream; //outputs eta to paramfile
  std::ofstream fstoutputstream;
  
  void checkLociNames(AdmixOptions *,InputData *data_);
  
  void getLabels( const std::string, Vector_i, std::string* );
  
  void loadAlleleStatesAndDistances(std::vector<std::string>* ChrmLabels,AdmixOptions *, 
				    InputData *data_,LogWriter *);
  
  void OpenFSTFile(AdmixOptions *options,LogWriter *Log); 
  
  // this can be moved when LociCorrSummary is moved
  static double strangExp( double );
 
  // we have four different methods to initalize allele freqs
  // 2nd and 3rd take i th locus as argument
  // other two loop over all composite loci
  // would be simpler to have one method
  // these methods are called by method LoadAlleleFreqs  
  void InitialiseAlleleFreqs(Matrix_d NewAlleleFreqs, int Populations);
  void InitialisePriorAlleleFreqs(Matrix_d New, int i, bool fixed);
  void InitialiseHistoricAlleleFreqs(Matrix_d New, int i);
  void SetDefaultAlleleFreqs(int Pops);

  // separate functions to set and to get allele probs
  void InitializeAlleleProbs();
  void SetAlleleProbs();
  // x is allele number
  double GetAlleleProbs( int x, int ancestry , int locus);

  void SamplePriorAlleleFreqs1D( Vector_d eta , int );
  void SamplePriorAlleleFreqsMultiDim( Vector_d eta , int);
  void SampleAlleleFreqs( int );

  // why is this function private and also public? 
  // void UpdateAlleleCounts(const std::vector<unsigned int>&, Vector_i );

  void UpdateAlleleCounts_HaploidData(const std::vector<unsigned int>&, int );
  // should have just one method to update allele freqs
  void UpdatePriorAlleleFreqs( int, const std::vector<Vector_d>& );

};

// functions required to update proportion vector Mu with adaptive rejection sampler
// likelihood, 1st and 2nd derivatives of log-likelihood
double fMu( Vector_d &, MatrixArray_i &, MatrixArray_d &, double );
double dfMu( Vector_d &, MatrixArray_i &, MatrixArray_d &, double );
double ddfMu( Vector_d &, MatrixArray_i &, MatrixArray_d &, double );


#endif
