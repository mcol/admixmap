// *-*-C++-*-*
#ifndef ALLELEFREQS_H
#define ALLELEFREQS_H 1
#include "InputData.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"
#include "LocusVisitor.h"
#include "AlleleFreqOutputter.h"

class AlleleFreqs{

public:
  AlleleFreqs(Genome *pLoci);
  ~AlleleFreqs();
  void Initialise(AdmixOptions *options, const Matrix_d& etaprior,LogWriter *Log,std::string *PopulationLabels);
  void Update(int iteration,int);

  //initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log);

  //outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  void OutputErgodicAvg( int iteration,AdmixOptions *options, std::ofstream *avgstream);
  //output samples of dispersion parameters (eta) to dispparamfile
  void OutputEta(int iteration, AdmixOptions *options, std::ofstream *LogFileStreamPtr);

  Genome *getLoci();
  CompositeLocus *getLocus(int);

  int GetNumberOfCompositeLoci();

  void OutputAlleleFreqs();

  void OutputFST(bool IsPedFile);

  void LoadAlleleFreqs(AdmixOptions *options, Chromosome ***chrm,LogWriter *Log, InputData *data,std::string **PopulationLabels);

  void ResetAlleleCounts();//resets Allelecounts to zero at start of iteration
  int IsRandom();//possibly should be bool?
  void UpdateFst();
  double *GetStatsForEta( int , int locus);
  double GetAlleleProbsMAP( int x, int ancestry , int locus);
  Vector_d GetPriorAlleleFreqs( int locus, int population );
  Vector_i GetAlleleCounts( int locus, int population );
  Vector_d getAlleleFreqsMAP( int locus, int population );
  Matrix_d &GetAlleleFreqs(int locus);
  Matrix_d *GetAlleleFreqs();
  Matrix_i &GetAlleleCounts(int locus);
  
  Matrix_d AlleleFreqs::GetSumAlleleFreqs(int locus);//is this used?

  void UpdateAlleleCounts(int locus, int h[2], Vector_i ancestry );
  void UpdateAlleleCounts_HaploidData(int locus, unsigned short **genotype, int ancestry );
  void ResetSumAlleleFreqs();
  void setAlleleFreqsMAP();
 
  void GetGenotypeProbs( double **Probs, int locus, unsigned short **genotype, 
			 std::vector<hapPair > &Haplotypes, bool diploid, bool fixed);
  void GetGenotypeProbs(double **Prob, int locus, unsigned short *genotype, Vector_i Haplotypes, bool diploid, bool fixed);

 // function to merge rare haplotypes for construction of score tests
  void SetMergedHaplotypes(Vector_d *alpha0, std::ofstream *LogFileStreamPtr, bool IsPedFile);

private:
  int Populations, NumberOfCompositeLoci;
  double *eta; //dispersion parameter
  double *SumEta;
  double *psi,*tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
 
  Matrix_d *Freqs;// allele frequencies except for last allele
  Matrix_d *AlleleFreqsMAP; // posterior mode of allele freqs
  Matrix_d *HistoricAlleleFreqs;
  Matrix_i *AlleleCounts;
  Matrix_d *HistoricLikelihoodAlleleFreqs;
  Matrix_d *PriorAlleleFreqs;

  Matrix_d *SumAlleleFreqs;// used to compute ergodic average

  dmatrix Fst;
  dmatrix SumFst;
  bool IsHistoricAlleleFreq;//indicator for dispersion model
  int RandomAlleleFreqs;//indicator for whether allele freqs are fixed or random - should be bool?

  Genome *Loci;//pointer to Loci object

  TuneRW *TuneEtaSampler;
  int Number,w; // Number is the number of updates of eta. The eta sampler is tuned every w updates. 

  double *etastep;
  double etastep0;

  int *NumberAccepted;
  double *SumAcceptanceProb;

  double *pp;//used to set merged haplotypes, which are used in the allelic association test

//    DARS SampleMu;
   std::vector<TuneRW> *MuProposal;

  LocusVisitor* allelefreqoutput;// object to output allele frequencies
  std::ofstream outputstream;//outputs eta to paramfile
  std::ofstream fstoutputstream;

  void checkLociNames(AdmixOptions *,InputData *data_);
    
  void getLabels( const std::string, Vector_i, std::string* );
 
  void OpenFSTFile(AdmixOptions *options,LogWriter *Log); 

  // we have four different functions to initialize allele freqs
  // 2nd and 3rd take i th locus as argument
  // other two loop over all composite loci
  // would be simpler to have one method
  // these functions are called by method LoadAlleleFreqs  
  void InitialiseAlleleFreqs(Matrix_d NewAlleleFreqs, int i, int Populations);
  void InitialisePriorAlleleFreqs(Matrix_d New, int i, bool fixed);
  void InitialiseHistoricAlleleFreqs(Matrix_d New, int i);
  void SetDefaultAlleleFreqs(int Pops);

  void SamplePriorAlleleFreqs1D( int );
  void SamplePriorAlleleFreqsMultiDim( int);
  void SampleAlleleFreqs(int, int);
  void UpdatePriorAlleleFreqs( int, const std::vector<Vector_d>& );
 
};
// functions required to update proportion vector Mu with adaptive rejection sampler
// likelihood, 1st and 2nd derivatives of log-likelihood
//Note that these are not part of AlleleFreqs class
double fMu( Vector_d &, Matrix_i &, Matrix_d &, double );
double dfMu( Vector_d &, Matrix_i &, Matrix_d &, double );
double ddfMu( Vector_d &, Matrix_i &, Matrix_d &, double );





#endif
