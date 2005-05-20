#include "IndividualCollection.h"
#include "StringSplitter.h"
#include "Regression.h"

using namespace std;

IndividualCollection::IndividualCollection()
{   
}

IndividualCollection::~IndividualCollection()
{
  //for(unsigned int i=0;i<_child.size();i++){
  for(unsigned int i=0;i<NumInd;i++){
    delete _child[i];
  }
  delete[] _child;
  delete indadmixoutput;
}

IndividualCollection::IndividualCollection(AdmixOptions* options,InputData *Data, Genome& Loci, Chromosome **chrm)
{
  Vector_i null_Vector_i(1);
  Matrix_d nullMatrix(1,1);
  OutcomeType = null_Vector_i;
  indadmixoutput = 0;
  TargetLabels = 0;
  LogLikelihood=0.0;
  SumLogLikelihood = 0.0;
  Covariates = nullMatrix;
  NumInd = Data->getNumberOfIndividuals();

  _child = new Individual*[NumInd];
  Individual::SetStaticMembers(Loci.GetNumberOfChromosomes(), &Loci);
    // Fill separate individuals.
  for (unsigned int i = 0; i < NumInd; ++i) {
    _child[i] = new Individual(i+1,options, Data, Loci, chrm);
    }
 
if (options->getAnalysisTypeIndicator() == 2)
    {
      OutcomeType.SetNumberOfElements( 1 );
      OutcomeType(0) = 0;
    }
  else if (options->getAnalysisTypeIndicator() == 3)
    {
      OutcomeType.SetNumberOfElements( 1 );
      OutcomeType(0) = 1;
    }
  else if (options->getAnalysisTypeIndicator() == 4)
    {
      OutcomeType.SetNumberOfElements( 1 );
      OutcomeType(0) = 1;
    }

}

void
IndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  for(unsigned int i=0; i<NumInd; i++){
    indadmixoutput->visitIndividual(*_child[i], _locusfortest, LogLikelihood);
  }
}

int
IndividualCollection::getSize()
{
  return NumInd;
}

Vector_i
IndividualCollection::GetSumXi()
{
   Vector_i sumxi = (*_child[0]).getSumXi();
   for( unsigned int i = 1; i <NumInd; i++ )
      sumxi += (*_child[i]).getSumXi();
   return sumxi;
}

double
IndividualCollection::GetSumrho0()
{
   double Sumrho0 = 0;
   for( unsigned int i = 0; i < NumInd; i++ )
      Sumrho0 += (*_child[i]).getSumrho0();
   return Sumrho0;
}

double
IndividualCollection::GetSumrho()
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < NumInd; i++ )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

MatrixArray_d IndividualCollection::getOutcome(){
  return Target;
}
Matrix_d IndividualCollection::getOutcome(int j){
  return Target(j);
}
Vector_d IndividualCollection::getTargetCol(int j, int k){
  return Target(j).GetColumn(k);
}
int IndividualCollection::getTargetSize(){
  return Target.GetNumberOfElements();
}

int IndividualCollection::getOutcomeType(int i){
  return OutcomeType(i);
}
Vector_i *IndividualCollection::getOutcomeType(){
  return &OutcomeType;
  }

Individual* IndividualCollection::getIndividual(int num)
{
  if (num < (int)NumInd){
    return _child[num];
  } else {
    return 0;
  }
}

void IndividualCollection::setAdmixtureProps(Matrix_d a)
{
  for(unsigned int i=0; i<NumInd; i++){
    _child[i]->setAdmixtureProps(a);
  }
}

void IndividualCollection::setAdmixturePropsX(Matrix_d a)
{
  for(unsigned int i=0; i<NumInd; i++){
    _child[i]->setAdmixturePropsX(a);
  }
}

MatrixArray_d IndividualCollection::getAncestries()
{
  MatrixArray_d ancestries(NumInd);
  for(unsigned int i=0; i<NumInd; i++){
    ancestries(i) = _child[i]->getAdmixtureProps();
  }
  return ancestries;
}

int IndividualCollection::GetNumberOfInputRows(){
  return Input.GetNumberOfRows();
}
int IndividualCollection::GetNumberOfInputCols(){
  return Input.GetNumberOfCols();
}
Matrix_d IndividualCollection::getCovariates(){
  return Covariates;
}

std::string IndividualCollection::getCovariateLabels(int i){
  return CovariateLabels[i];
  }
std::string *IndividualCollection::getCovariateLabels(){
  return CovariateLabels;
  }

double IndividualCollection::getExpectedY(int i){
  return ExpectedY(0)(i,0);
 }

std::string IndividualCollection::getTargetLabels(int k){
  return TargetLabels[k];
}

void IndividualCollection::SetExpectedY(int k,Matrix_d beta){
  ExpectedY(k) = Covariates * beta;
}
void IndividualCollection::calculateExpectedY( int k)
{
  for(unsigned int i = 0; i < NumInd; i++ )
    ExpectedY(k)( i, 0 ) = 1 / ( 1 + exp( -ExpectedY(k)( i, 0 ) ) );
}

void IndividualCollection::Initialise(AdmixOptions *options,MatrixArray_d *beta, Genome *Loci, std::string *PopulationLabels, 
				      double rhoalpha, double rhobeta, LogWriter *Log, const Matrix_d &MLEMatrix){
  //Open indadmixture file  
  if ( strlen( options->getIndAdmixtureFilename() ) ){
    Log->logmsg(true,"Writing individual-level parameters to ");
    Log->logmsg(true,options->getIndAdmixtureFilename());
    Log->logmsg(true,"\n");
    indadmixoutput = new IndAdmixOutputter(options,Loci,PopulationLabels);
  }
  else {
    Log->logmsg(true,"No indadmixturefile given\n");
  }
  //Set locusfortest if specified
 if( options->getLocusForTestIndicator() )
     _locusfortest = Loci->GetChrmAndLocus( options->getLocusForTest() );

 //Initialise Admixture Proportions
  Matrix_d admix_null;
   if( options->isRandomMatingModel() )
    admix_null.SetNumberOfElements(options->getPopulations(),2);
  else
    admix_null.SetNumberOfElements(options->getPopulations(),1);
  admix_null.SetElements( (double)1.0 / options->getPopulations() );
  Vector_d alphatemp;

  if( options->sizeInitAlpha() == 0 ){
    alphatemp.SetNumberOfElements( options->getPopulations() );
    alphatemp.SetElements( 1.0 );
  }
  else if( options->sizeInitAlpha() == 1 ){
    alphatemp = options->getInitAlpha(0);
  }
  else if( options->getAnalysisTypeIndicator() < 0 ){
    alphatemp = options->getInitAlpha(0);
 
    for( int k = 0; k < options->getPopulations(); k++ ){
       if( alphatemp(k) == 0 ) admix_null(k,0) = 0.0;
     }
     
     alphatemp = options->getInitAlpha(1);
  
     for( int k = 0; k < options->getPopulations(); k++ ){
       if( alphatemp(k) == 0 ) admix_null(k,1) = 0.0;
     }
  }
  setAdmixtureProps(admix_null);
  if( Loci->isX_data() )setAdmixturePropsX(admix_null);

  //Regression stuff
  if(options->getAnalysisTypeIndicator() >=2){
    ExpectedY.SetNumberOfElementsWithDimensions( getTargetSize(), NumInd, 1 );
    //Covariates.SetNumberOfElements(1);
    Matrix_d temporary( NumInd, 1 );
    temporary.SetElements(1);
      
    if( Input.GetNumberOfRows() == (int)NumInd ){
      Covariates = ConcatenateHorizontally( temporary, Input );
      Vector_d mean;
      mean = Input.ColumnMean();
      for(unsigned int i = 0; i < NumInd; i++ )
	for( int j = 0; j < Input.GetNumberOfCols(); j++ )
	  Input( i, j ) -= mean(j);
    } else {
      Covariates = temporary;
    }

    
    if( !options->getScoreTestIndicator() && options->getPopulations() > 1 ){
      temporary.SetNumberOfElements( NumInd, options->getPopulations() - 1 );
      temporary.SetElements( 1 / options->getPopulations() );
      Covariates = ConcatenateHorizontally( Covariates, temporary );
    }

    for( int k = 0; k < getTargetSize(); k++ ){
      SetExpectedY(k,(*beta)(k));
      if( getOutcomeType(k) )calculateExpectedY(k);
    }
  }
  //Misc.
  SumLogTheta.SetNumberOfElements( options->getPopulations());
  InitialiseMLEs(rhoalpha,rhobeta,options, MLEMatrix);
  //set to very large negative value (effectively -Inf) so the first value is guaranteed to be greater
  MaxLogLikelihood.assign(NumInd, -9999999 );
}

void IndividualCollection::InitialiseMLEs(double rhoalpha, double rhobeta, AdmixOptions * options, const Matrix_d &MLEMatrix){
  //set thetahat and rhohat, estimates of individual admixture and sumintensities
   thetahat.SetNumberOfElementsWithDimensions( NumInd, 1, 1 );
   thetahatX.SetNumberOfElementsWithDimensions( NumInd, 1, 1 );

   vector<double> r(2, rhoalpha/rhobeta );
   rhohat.resize( NumInd, r );
   rhohatX.resize( NumInd, r );

   //use previously read values from file, if available
   if( options->getAnalysisTypeIndicator() == -2 ){
      rhohat[0][0] = MLEMatrix( options->getPopulations(), 0 );
      if( options->getXOnlyAnalysis() )
	thetahat(0) = MLEMatrix.SubMatrix( 0, options->getPopulations() - 1, 0, 0 );
      else{
	thetahat(0) = MLEMatrix.SubMatrix( 0, options->getPopulations() - 1, 0, 1 );
	rhohat[0][1] = MLEMatrix(options->getPopulations(), 1 );
      }
      setAdmixtureProps(thetahat(0));
   }
   else if( options->getAnalysisTypeIndicator() == -1 ){
      thetahat = getAncestries();
   }
 
}

void IndividualCollection::LoadGenotypes(AdmixOptions *options, InputData *data_, LogWriter *Log, Genome *Loci){
  // Load Genotypes into Individuals

  /**This block should be replaced or removed. Omitted by Dmitry
  if(options->IsPedFile() == 0){
    Matrix_d locifileData;
    locifileData.Load( options->getGeneInfoFilename() );
    lshtm_match m;
    lshtm_path_regex("\"[[:space:]]+\"", &m, options->getGeneticDataFilename(), REG_EXTENDED);
    for( unsigned int i = 0; i < m.lCount; i++ ){
      if( (int)(m.mPerLine[i]+1-options->genotypesSexColumn()) != locifileData.GetNumberOfRows() && m.mPerLine[i] != 0 ){
	Log->logmsg(true,"Error: Line ");
	Log->logmsg(true, (int)i+1);
	Log->logmsg(true," of the genotypes file ");
	Log->logmsg(true, options->getGeneticDataFilename());
	Log->logmsg(true," has ");
	Log->logmsg(true,(int)(m.mPerLine[i]+1-options->genotypesSexColumn()));
	Log->logmsg(true," columns\n");
	Log->logmsg(true,"But the loci file ");
	Log->logmsg(true,options->getGeneInfoFilename());
	Log->logmsg(true," has ");
	Log->logmsg(true,locifileData.GetNumberOfRows());
	Log->logmsg(true," rows\n");
	exit(0);
      }
    }
  }
  */  
  if ( strlen( options->getInputFilename() ) != 0 ){  
    LoadCovariates(options,data_, Log);
  }
  if ( Input.GetNumberOfMissingValues() ) Input.SetMissingValuesToColumnMeans();

  if ( strlen( options->getTargetFilename() ) != 0 ){
    LoadOutcomeVar(options, data_, Log);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(options, data_, Log);
  }
  CheckGenotypes(Loci, Log);
}

void IndividualCollection::LoadCovariates(AdmixOptions *options, InputData *data_, LogWriter *Log){
  //LOAD INPUT (COVARIATES)

  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getInputFilename());
  Log->logmsg(false,".\n");
  Input = data_->getInputMatrix();

  if( (int)NumInd != Input.GetNumberOfRows() - 1 ){
    Log->logmsg(true,"Genotype file has ");
    Log->logmsg(true,NumInd);
    Log->logmsg(true," observations and Input file has ");
    Log->logmsg(true,Input.GetNumberOfRows() - 1);
    Log->logmsg(true," observations.\n");
    exit(0);
  }
  if( options->getTextIndicator() ){
    Input.SubMatrix2( 1, NumInd, 0, Input.GetNumberOfCols() - 1 );
    CovariateLabels = new string[ Input.GetNumberOfCols() ];
    Vector_i vtemp( Input.GetNumberOfCols() );
    vtemp.SetElements( 1 );
    getLabels(data_->getInputData()[0], vtemp, CovariateLabels);//?getCovariatelabels()
  }
}

void IndividualCollection::LoadOutcomeVar(AdmixOptions *options, InputData *data_, LogWriter *Log){
  //LOAD TARGET (OUTCOME VARIABLE)
  Matrix_d TempTarget, temporary;
  string *TempLabels = 0;

  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getTargetFilename());
  Log->logmsg(false,".\n");
  //const Matrix_d& LoadTarget = data_->getTargetMatrix();
  //conversion necessary because LoadTarget is changed further down
  Matrix_d& LoadTarget = (Matrix_d&)data_->getTargetMatrix();
  if( LoadTarget.GetNumberOfRows() - 1 != (int)NumInd ){
    Log->logmsg(true,"Inconsistency in number of rows in outcomevarfile and genotypefile.\n");
  }
  TempLabels = new string[ LoadTarget.GetNumberOfCols() ];

  Vector_i vtemp( LoadTarget.GetNumberOfCols() );
  vtemp.SetElements(1);
  getLabels(data_->getTargetData()[0], vtemp, TempLabels);

  if( options->getAnalysisTypeIndicator() == 5 ){
    TargetLabels = new string[ LoadTarget.GetNumberOfCols() ];
    Target.SetNumberOfElements(LoadTarget.GetNumberOfCols());
    OutcomeType.SetNumberOfElements( LoadTarget.GetNumberOfCols() );

    for( int j = 0; j < LoadTarget.GetNumberOfCols(); j++ ){
      TargetLabels[j] = TempLabels[j];
      TempTarget = LoadTarget;
      TempTarget.SubMatrix2( 1, NumInd, j, j );
      for(unsigned int i = 0; i < NumInd; i++ ){
	if( !TempTarget.IsMissingValue( i, 0 ) &&
	    (TempTarget( i, 0 ) == 0 || TempTarget( i, 0 ) == 1) )
	  OutcomeType(j) = 1;
      }
      Target(j) = TempTarget;

      if( getOutcomeType(j) )
	{
	  Log->logmsg(true,"Binary variable: ");
	  Log->logmsg(true,getTargetLabels(j));
	  Log->logmsg(true,".\n");
	}
      else
	{
	  Log->logmsg(true,"Continuous variable: ");
	  Log->logmsg(true,getTargetLabels(j));
	  Log->logmsg(true,".\n");
	}
    }
  }
  else{
    TargetLabels = new string[ 1 ];
    TargetLabels[0] = TempLabels[ options->getTargetIndicator() ];
    Target.SetNumberOfElements(1);
    Log->logmsg(true,"Regressing on: ");
    Log->logmsg(true,getTargetLabels(0));
    Log->logmsg(true,".\n");
    if( LoadTarget.GetNumberOfRows() - 1 != (int)NumInd ){
      Log->logmsg(true,"Target file has ");
      Log->logmsg(true,LoadTarget.GetNumberOfRows() - 1);
      Log->logmsg(true," observations and Genotypes file has ");
      Log->logmsg(true,NumInd);
      Log->logmsg(true," observations.\n");
      exit(1);
    }
    LoadTarget.SubMatrix2( 1, NumInd, options->getTargetIndicator(), options->getTargetIndicator() );
    Target(0) = LoadTarget;
  }

  delete [] TempLabels;
}

void IndividualCollection::LoadRepAncestry(AdmixOptions *options, InputData *data_, LogWriter *Log){
  //LOAD REPORTED ANCESTRY IF GIVEN   

  ReportedAncestry.SetNumberOfElements( NumInd );
  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getReportedAncestryFilename());
  Log->logmsg(false,".\n");
  const Matrix_d& temporary = data_->getReportedAncestryMatrix();

  if( temporary.GetNumberOfRows() != 2 * (int)NumInd ){
    Log->logmsg(false,"Error\n");
    Log->logmsg(false,options->getReportedAncestryFilename());
    Log->logmsg(false," has ");
    Log->logmsg(false,temporary.GetNumberOfRows());
    Log->logmsg(false," rows\n");
    Log->logmsg(false,options->getGeneticDataFilename());
    Log->logmsg(false," has ");
    Log->logmsg(false,NumInd);
    Log->logmsg(false," rows\n");
    exit(0);}
  if( temporary.GetNumberOfCols() != options->getPopulations() ){
    Log->logmsg(false,"Error\n");
    Log->logmsg(false,options->getReportedAncestryFilename());
    Log->logmsg(false," has ");
    Log->logmsg(false,temporary.GetNumberOfCols());
    Log->logmsg(false," cols\n");
    Log->logmsg(false,options->getAlleleFreqFilename());
    Log->logmsg(false," has ");
    Log->logmsg(false,options->getPopulations());
    Log->logmsg(false," cols\n");
    exit(0);
  }
  for( int i = 0; i < temporary.GetNumberOfRows() / 2; i++ )
    ReportedAncestry(i) = temporary.SubMatrix( 2*i, 2*i + 1, 0, temporary.GetNumberOfCols() - 1 );
 
}

void IndividualCollection::getLabels( const string buffer, Vector_i temporary, string *labels )
{
  StringSplitter splitter;
  const Vector_s& labels_tmp = splitter.split(buffer);

  for (size_t i = 0, index = 0; i < labels_tmp.size(); ++i) {
    if (temporary.GetNumberOfElements() == 1 || temporary(i)) {            
      labels[index++] = labels_tmp[i];
    }
  }
}
void IndividualCollection::getLabels(const Vector_s& data, Vector_i temporary, string *labels)
{
    for (size_t i = 0, index = 0; i < data.size(); ++i) {
        if (temporary.GetNumberOfElements() == 1 || temporary(i)) {            
            labels[index++] = data[i];
        }
    }
}

//should be in InputData class
void IndividualCollection::CheckGenotypes(Genome *Loci,LogWriter *Log)
{
  bool error = false;
   
  for(unsigned int i = 0; i < NumInd; i++ ){
    for(unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      for( int k = 0; k < (*Loci)(j)->GetNumberOfLoci(); k++ ){
	Individual* ind =  _child[i];
	if( (int)(ind->getGenotype(j)[k][0])   > (*Loci)(j)->GetNumberOfAllelesOfLocus( k ) || 
	    (int)(ind->getGenotype(j)[k][1]) > (*Loci)(j)->GetNumberOfAllelesOfLocus( k ) ){
	  Log->logmsg(false, "Error in genotypes file\n");
	  Log->logmsg(false, "Individual ");
	  Log->logmsg(false, i);
	  Log->logmsg(false, " at locus ");
	  Log->logmsg(false, (*Loci)(j)->GetLabel(k));
	  Log->logmsg(false, " has genotype ");
	  Log->logmsg(false,ind->getGenotype(j)[k][0]);Log->logmsg(false, " ");
	  Log->logmsg(false, ind->getGenotype(j)[k][1]);Log->logmsg(false, " \n");
	  Log->logmsg(false, "Number of allelic states at locus = ");
	  Log->logmsg(false, (*Loci)(j)->GetNumberOfAllelesOfLocus( k ));Log->logmsg(false, "\n");
	  error = true;
	}
      }
    }
  }
  if( error )
    exit(0);
}

void IndividualCollection::Update(int iteration, AlleleFreqs *A, Regression *R, Vector_d &poptheta, AdmixOptions *options,
				  Chromosome **chrm, vector<Vector_d> alpha, bool _symmetric, 
				  vector<bool> _admixed, double rhoalpha, double rhobeta,
				  std::ofstream *LogFileStreamPtr, chib *MargLikelihood){
  SumLogTheta.SetElements( 0.0 );
  if(iteration > options->getBurnIn())Individual::ResetScores(options);


  for(unsigned int i = 0; i < NumInd; i++ ){
    
    if( options->getPopulations() > 1 ){
      _child[i]->SampleParameters(i, &SumLogTheta, A, iteration , &Target, OutcomeType, ExpectedY, *(R->getlambda()), 
				  R->getNoCovariates(),  Covariates,*(R->getbeta()),poptheta, options, 
				  chrm, alpha, _symmetric, _admixed, rhoalpha, rhobeta, sigma,  
				  DerivativeInverseLinkFunction(options->getAnalysisTypeIndicator(), i),
				  R->getDispersion(OutcomeType(0)));}
    
    else{
      _child[i]->OnePopulationUpdate(i, &Target, OutcomeType, ExpectedY, *(R->getlambda()), options->getAnalysisTypeIndicator());
    }   
    
    if( options->getAnalysisTypeIndicator() < 0 && i == 0 )//check if this condition is correct
      _child[i]->ChibLikelihood(i,iteration, &LogLikelihood, &SumLogLikelihood, MaxLogLikelihood, 
				options, chrm, alpha,_admixed, rhoalpha, rhobeta,
				thetahat, thetahatX, rhohat, rhohatX,LogFileStreamPtr, MargLikelihood, A);
  }

}

void IndividualCollection::Output(std::ofstream *LogFileStreamPtr){
  //Used only for IndAdmixHierModel = 0
  for(unsigned  int i = 0; i < NumInd; i++ )
     *LogFileStreamPtr << thetahat(i).GetColumn(0) << " "
     << thetahat(i).GetColumn(1) << " "
     << rhohat[i][0] << " " << rhohat[i][1] << endl;
}

void IndividualCollection::OutputErgodicAvg(int samples, chib *MargLikelihood, std::ofstream *avgstream){
     *avgstream << SumLogLikelihood / samples << " "
               << MargLikelihood->getLogPosterior();
}

void
IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter *Log, AlleleFreqs *A, std::string *PopulationLabels)
{
   Log->logmsg(true,"Log-likelihood for unadmixed ");
   Log->logmsg(true, (*PopulationLabels)[0]);
   Log->logmsg(true, ": ");
   Log->logmsg(true, _child[0]->getLogLikelihoodOnePop(A));
   Log->logmsg(true,"\n");

}
double IndividualCollection::getSumLogTheta(int i){
  return SumLogTheta(i);
}
double IndividualCollection::getLL(){
  //not currently used
  return SumLogLikelihood;
}

//returns Derivative of Inverse Link Function for individual i
double IndividualCollection::DerivativeInverseLinkFunction(int AnalysisType,int i){
  double DInvLink = 1.0;
  double EY = getExpectedY(i);
  int OutcomeType = getOutcomeType(0);

    //Linear regression
    if(AnalysisType == 2 ){
      DInvLink = 1.0;
      }
    //Logistic Regression
    else if( AnalysisType == 3 || AnalysisType == 4 ){
      DInvLink = EY * (1.0 - EY);
    }
    else if( AnalysisType == 5 ){
      DInvLink = OutcomeType ? EY*(1.0-EY):1.0;
    }
 
  return DInvLink;    
}
