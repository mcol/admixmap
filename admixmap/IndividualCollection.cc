#include "IndividualCollection.h"
#include "StringSplitter.h"

using namespace std;

IndividualCollection::IndividualCollection()
{   
}

IndividualCollection::~IndividualCollection()
{
  for(unsigned int i=0;i<_child.size();i++){
    delete _child[i];
  }
  delete indadmixoutput;
}

IndividualCollection::IndividualCollection(AdmixOptions* options,
    const Matrix_s& data, Genome& Loci, Genome& chrm)
{
  Vector_i null_Vector_i(1);
  TargetType = null_Vector_i;
  indadmixoutput = 0;
  TargetLabels = 0;
  LogLikelihood=0.0;
  SumLogLikelihood = 0.0;

    // Determine if ped file.
    const int isPedFile = 2*data[0].size() - 1 == data[1].size() ? 1 : 0;
    options->IsPedFile(isPedFile);

    // Fill separate individuals.
    for (size_t i = 1; i < data.size() ; ++i) {
        _child.push_back(new Individual(options, data[i], Loci, chrm));
    }

if (options->getAnalysisTypeIndicator() == 2)
    {
      SetNumberOfTargetTypeElements(1);//TargetType->SetNumberOfElements( 1 );
      SetTargetType(0,0);//(*TargetType)( 0 ) = 0;
    }
  else if (options->getAnalysisTypeIndicator() == 3)
    {
      SetNumberOfTargetTypeElements(1);//TargetType->SetNumberOfElements( 1 );
      SetTargetType(0,1);//(*TargetType)( 0 ) = 1;
    }
  else if (options->getAnalysisTypeIndicator() == 4)
    {
      SetNumberOfTargetTypeElements(1);//TargetType->SetNumberOfElements( 1 );
      SetTargetType(0,1);//(*TargetType)( 0 ) = 1;
    }

}

void
IndividualCollection::accept()
{
  indadmixoutput->visitIndividualCollection(*this);
  for(unsigned int i=0; i<_child.size(); i++){
    _child[i]->accept( *indadmixoutput, ExpectedY(0)(i,0), _locusfortest, LogLikelihood );
  }
}

int
IndividualCollection::getSize()
{
  return _child.size();
}

Vector_i
IndividualCollection::GetSumXi()
{
   Vector_i sumxi = (*_child[0]).getSumXi();
   for( unsigned int i = 1; i < _child.size(); i++ )
      sumxi += (*_child[i]).getSumXi();
   return sumxi;
}

double
IndividualCollection::GetSumrho0()
{
   double Sumrho0 = 0;
   for( unsigned int i = 0; i < _child.size(); i++ )
      Sumrho0 += (*_child[i]).getSumrho0();
   return Sumrho0;
}

double
IndividualCollection::GetSumrho()
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < _child.size(); i++ )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

MatrixArray_d IndividualCollection::getTarget(){
  return Target;
}
Matrix_d IndividualCollection::getTarget(int j){
  return Target(j);
}
Vector_d IndividualCollection::getTargetCol(int j, int k){
  return Target(j).GetColumn(k);
}
int IndividualCollection::getTargetSize(){
  return Target.GetNumberOfElements();
}

int IndividualCollection::getTargetType(int i){
  return TargetType(i);
}
Vector_i *IndividualCollection::getTargetType(){
  return &TargetType;
  }

void
IndividualCollection::add(Individual* ind)
{
  _child.push_back(ind);
}

Individual*
IndividualCollection::getIndividual(int num)
{
  if (num < (int)_child.size()){
    return _child[num];
  } else {
    return 0;
  }
}

void
IndividualCollection::setAncestry(Matrix_d ancestry)
{
  for(unsigned int i=0; i<_child.size(); i++){
    _child[i]->setAncestry(ancestry);
  }
}

void
IndividualCollection::setAncestryX(Matrix_d ancestry)
{
  for(unsigned int i=0; i<_child.size(); i++){
    _child[i]->setAncestryX(ancestry);
  }
}

MatrixArray_d
IndividualCollection::getAncestries()
{
  MatrixArray_d ancestries(_child.size());
  for(unsigned int i=0; i<_child.size(); i++){
    ancestries(i) = _child[i]->getAncestry();
  }
  return ancestries;
}

int IndividualCollection::GetNumberOfInputRows(){
  return Input.GetNumberOfRows();
}
int IndividualCollection::GetNumberOfInputCols(){
  return Input.GetNumberOfCols();
}
MatrixArray_d IndividualCollection::getCovariates(){
  return Covariates;
}
Matrix_d IndividualCollection::getCovariates(int i){
  return Covariates(i);
}
std::string IndividualCollection::getCovariateLabels(int i){
  return CovariateLabels[i];
  }
std::string *IndividualCollection::getCovariateLabels(){
  return CovariateLabels;
  }

Matrix_d* IndividualCollection::getExpectedY0(){
   return &ExpectedY(0);
 }

void IndividualCollection::SetNumberOfTargetTypeElements(int i){
  TargetType.SetNumberOfElements(i);
}
void IndividualCollection::SetTargetType(int index, int setto){
  TargetType(index) = setto;
}

std::string IndividualCollection::getTargetLabels(int k){
  return TargetLabels[k];
}

void IndividualCollection::SetExpectedY(int k,Matrix_d beta){
ExpectedY(k) = getCovariates(0) * beta;
}
void IndividualCollection::calculateExpectedY( int k)
{
  for( int i = 0; i < getSize(); i++ )
    ExpectedY(k)( i, 0 ) = 1 / ( 1 + exp( -ExpectedY(k)( i, 0 ) ) );
}

void IndividualCollection::Initialise(AdmixOptions *options,MatrixArray_d *beta, Genome *Loci, std::string *PopulationLabels){

if ( strlen( options->getIndAdmixtureFilename() ) ){
    indadmixoutput = new IndAdmixOutputter(options,Loci,PopulationLabels);
  }
 if( options->getLocusForTestIndicator() )
     _locusfortest = Loci->GetChrmAndLocus( options->getLocusForTest() );

  Matrix_d ancestry_null;
   if( options->getModelIndicator() )
    ancestry_null.SetNumberOfElements(options->getPopulations(),2);
  else
    ancestry_null.SetNumberOfElements(options->getPopulations(),1);
  ancestry_null.SetElements( (double)1.0 / options->getPopulations() );
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
       if( alphatemp(k) == 0 ) ancestry_null(k,0) = 0.0;
     }
     
     alphatemp = options->getInitAlpha(1);
  
     for( int k = 0; k < options->getPopulations(); k++ ){
       if( alphatemp(k) == 0 ) ancestry_null(k,1) = 0.0;
     }
  }
  setAncestry(ancestry_null);
  if( Loci->isX_data() ){
    setAncestryX(ancestry_null);
  }

  if(options->getAnalysisTypeIndicator() >=2){
    ExpectedY.SetNumberOfElementsWithDimensions( getTargetSize(), getSize(), 1 );
    Covariates.SetNumberOfElements(1);
    Matrix_d temporary( getSize(), 1 );
    temporary.SetElements(1);
      
    if( Input.GetNumberOfRows() == getSize() ){
      Covariates(0) = ConcatenateHorizontally( temporary, Input );
      Vector_d mean;
      mean = Input.ColumnMean();
      for( int i = 0; i < getSize(); i++ )
	for( int j = 0; j < Input.GetNumberOfCols(); j++ )
	  Input( i, j ) -= mean(j);
    } else {
      Covariates(0) = temporary;
    }

    
    if( !options->getScoreTestIndicator() && options->getPopulations() > 1 ){
      temporary.SetNumberOfElements( getSize(), options->getPopulations() - 1 );
      temporary.SetElements( 1 / options->getPopulations() );
      Covariates(0) = ConcatenateHorizontally( Covariates(0), temporary );
    }

    for( int k = 0; k < getTargetSize(); k++ ){
      SetExpectedY(k,(*beta)(k));
      if( getTargetType(k) )
	calculateExpectedY(k);
    }
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

  if( getSize() != Input.GetNumberOfRows() - 1 ){
    Log->logmsg(true,"Genotype file has ");
    Log->logmsg(true,getSize());
    Log->logmsg(true," observations and Input file has ");
    Log->logmsg(true,Input.GetNumberOfRows() - 1);
    Log->logmsg(true," observations.\n");
    exit(0);
  }
  if( options->getTextIndicator() ){
    Input.SubMatrix2( 1, getSize(), 0, Input.GetNumberOfCols() - 1 );
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
  if( LoadTarget.GetNumberOfRows() - 1 != getSize() ){
    Log->logmsg(true,"Inconsistency in number of rows in targetfile and genotypefile.\n");
  }
  TempLabels = new string[ LoadTarget.GetNumberOfCols() ];

  Vector_i vtemp( LoadTarget.GetNumberOfCols() );
  vtemp.SetElements(1);
  getLabels(data_->getTargetData()[0], vtemp, TempLabels);

  if( options->getAnalysisTypeIndicator() == 5 ){
    TargetLabels = new string[ LoadTarget.GetNumberOfCols() ];
    Target.SetNumberOfElements(LoadTarget.GetNumberOfCols());
    TargetType.SetNumberOfElements( LoadTarget.GetNumberOfCols() );

    for( int j = 0; j < LoadTarget.GetNumberOfCols(); j++ ){
      TargetLabels[j] = TempLabels[j];
      TempTarget = LoadTarget;
      TempTarget.SubMatrix2( 1, getSize(), j, j );
      for( int i = 0; i < getSize(); i++ ){
	if( !TempTarget.IsMissingValue( i, 0 ) &&
	    (TempTarget( i, 0 ) == 0 || TempTarget( i, 0 ) == 1) )
	  TargetType(j) = 1;
      }
      Target(j) = TempTarget;

      if( getTargetType(j) )
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
    if( LoadTarget.GetNumberOfRows() - 1 != getSize() ){
      Log->logmsg(true,"Target file has ");
      Log->logmsg(true,LoadTarget.GetNumberOfRows() - 1);
      Log->logmsg(true," observations and Genotypes file has ");
      Log->logmsg(true,getSize());
      Log->logmsg(true," observations.\n");
      exit(1);
    }
    LoadTarget.SubMatrix2( 1, getSize(), options->getTargetIndicator(), options->getTargetIndicator() );
    Target(0) = LoadTarget;
  }

  delete [] TempLabels;
}

void IndividualCollection::LoadRepAncestry(AdmixOptions *options, InputData *data_, LogWriter *Log){
  //LOAD REPORTED ANCESTRY IF GIVEN   

  ReportedAncestry.SetNumberOfElements( getSize() );
  Log->logmsg(false,"Loading ");
  Log->logmsg(false,options->getReportedAncestryFilename());
  Log->logmsg(false,".\n");
  const Matrix_d& temporary = data_->getReportedAncestryMatrix();

  if( temporary.GetNumberOfRows() != 2 * getSize() ){
    Log->logmsg(false,"Error\n");
    Log->logmsg(false,options->getReportedAncestryFilename());
    Log->logmsg(false," has ");
    Log->logmsg(false,temporary.GetNumberOfRows());
    Log->logmsg(false," rows\n");
    Log->logmsg(false,options->getGeneticDataFilename());
    Log->logmsg(false," has ");
    Log->logmsg(false,getSize());
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

void IndividualCollection::CheckGenotypes(Genome *Loci,LogWriter *Log)
{
  bool error = false;
   
  for( int i = 0; i < getSize(); i++ ){
    for( int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      for( int k = 0; k < (*Loci)(j)->GetNumberOfLoci(); k++ ){
	Individual* ind = getIndividual(i);
	if( (int)(ind->getGenotype(j)[2*k])   > (*Loci)(j)->GetNumberOfAllelesOfLocus( k ) || 
	    (int)(ind->getGenotype(j)[2*k+1]) > (*Loci)(j)->GetNumberOfAllelesOfLocus( k ) ){
	  Log->logmsg(false, "Error in genotypes file\n");
	  Log->logmsg(false, "Individual ");
	  Log->logmsg(false, i);
	  Log->logmsg(false, " at locus ");
	  Log->logmsg(false, (*Loci)(j)->GetLabel(k));
	  Log->logmsg(false, " has genotype ");
	  Log->logmsg(false,ind->getGenotype(j)[2*k]);Log->logmsg(false, " ");
	  Log->logmsg(false, ind->getGenotype(j)[2*k+1]);Log->logmsg(false, " \n");
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

void IndividualCollection::Update(int iteration, Vector_d *SumLogTheta, Vector_d *lambda, int NoCovariates,  
				  MatrixArray_d *beta, Vector_d &poptheta, AdmixOptions *options,
				  Vector_d &f, Genome *Loci, Genome *chrm, vector<Vector_d> alpha, bool _symmetric, 
				  vector<bool> _admixed, double rhoalpha, double rhobeta,
				  std::ofstream *LogFileStreamPtr, chib *MargLikelihood){

for( int i = 0; i < getSize(); i++ ){
  getIndividual(i)->IndivUpdate(i,iteration, 
				SumLogTheta, &Target, TargetType, ExpectedY, *lambda,
				NoCovariates, Covariates(0), *beta, poptheta, options,
				f, Loci, chrm, alpha, _symmetric, 
				_admixed, rhoalpha, rhobeta, sigma);

  if( options->getAnalysisTypeIndicator() == -3 || 
      options->getAnalysisTypeIndicator() == -1 && i == 0 /*&& options->LikOutput()*/)
    getIndividual(i)->ChibLikelihood(i,iteration, &LogLikelihood, &SumLogLikelihood, MaxLogLikelihood, 
				     options, Loci, chrm, alpha,_admixed, rhoalpha, rhobeta,
				     thetahat, thetahatX, rhohat, rhohatX,LogFileStreamPtr, MargLikelihood);
 }
}

void IndividualCollection::PreUpdate(double rhoalpha, double rhobeta, AdmixOptions * options){
   Matrix_d temp;
   thetahat.SetNumberOfElementsWithDimensions( getSize(), 1, 1 );
   thetahatX.SetNumberOfElementsWithDimensions( getSize(), 1, 1 );

   vector<double> r(2, rhoalpha/rhobeta );
   rhohat.resize( getSize(), r );
   rhohatX.resize( getSize(), r );

   if( options->getAnalysisTypeIndicator() == -2 ){
      temp.Load( options->getMLEFilename() );
      rhohat[0][0] = temp( options->getPopulations(), 0 );
      if( options->getXOnlyAnalysis() )
	thetahat(0) = temp.SubMatrix( 0, options->getPopulations() - 1, 0, 0 );
      else{
	thetahat(0) = temp.SubMatrix( 0, options->getPopulations() - 1, 0, 1 );
	rhohat[0][1] = temp(options->getPopulations(), 1 );
      }
      setAncestry(thetahat(0));
   }
   else if( options->getAnalysisTypeIndicator() == -1 ){
      thetahat = getAncestries();
   }
 MaxLogLikelihood.assign(getSize(), -9999999 );
}

void IndividualCollection::Output(std::ofstream *LogFileStreamPtr){
  //Used only for AnalysisTypeIndicator = -3
  for( int i = 0; i < getSize(); i++ )
     *LogFileStreamPtr << thetahat(i).GetColumn(0) << " "
     << thetahat(i).GetColumn(1) << " "
     << rhohat[i][0] << " " << rhohat[i][1] << endl;
}

void IndividualCollection::OutputErgodicAvg(int samples, chib *MargLikelihood, std::ofstream *avgstream){
     *avgstream << SumLogLikelihood / samples << " "
               << MargLikelihood->getLogPosterior();
}

void
IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter *Log, Genome *Loci, std::string *PopulationLabels)
{
   Individual* ind = getIndividual(0);
   Log->logmsg(true,"Log-likelihood for unadmixed ");
   Log->logmsg(true, (*PopulationLabels)[0]);
   Log->logmsg(true, ": ");
   Log->logmsg(true,ind->getLogLikelihoodOnePop(*Loci));
   Log->logmsg(true,"\n");

}
double IndividualCollection::getLL(){
  return SumLogLikelihood;
}
