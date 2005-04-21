#include "Individual.h"
#include "StringConvertor.h"

#define PR(x) cout << #x << " = " << x << endl;

Matrix_d Individual::AffectedsScore;
Matrix_d Individual::AffectedsVarScore;
Matrix_d Individual::AffectedsInfo;
MatrixArray_d Individual::AncestryScore;
MatrixArray_d Individual::AncestryInfo;
Matrix_d Individual::AncestryVarScore;
Matrix_d Individual::AncestryInfoCorrection;
Matrix_d Individual::B;

Individual::Individual()
{
}

Individual::Individual(int mynumber,AdmixOptions* options, InputData *Data, Genome& Loci,Chromosome **chrm)
{

    if( options->getRhoIndicator() ){
        TruncationPt = options->getTruncPt();
        if( options->getModelIndicator() )
            if(options->getRho() < 90.0 )
                _rho.assign(2,options->getRho());
            else
                _rho.assign(2,1);
        else
            if(options->getRho() < 90.0 )
                _rho.assign(1,options->getRho());
            else
                _rho.assign(2,1);
    }
  for(int j=0 ;j<2;++j) f[j] = new double[Loci.GetNumberOfCompositeLoci()];

    numChromosomes = Loci.GetNumberOfChromosomes();
    //numCompLoci = Loci.GetSizesOfChromosomes();

    // Read sex value if present.
    if (options->genotypesSexColumn() == 1) {
	sex = Data->GetSexValue(mynumber);
    }

    int numCompositeLoci = Loci.GetNumberOfCompositeLoci();

    //create vector of jump indicators
    vector<bool> true_vector(numCompositeLoci,true);
    _xi.assign(2,true_vector);
    sumxi.SetNumberOfElements(numCompositeLoci);

    PossibleHaplotypes = new Vector_i[numCompositeLoci]; // vector of possible haplotype pairs - expect 2 integers per locus 
    LocusAncestry.SetNumberOfElements( Loci.GetNumberOfChromosomes() ); // array of matrices in which each col stores 2 integers 
 
    AncestryProbs = 0;                                                                      // or 1 integer (haploid) 
    if( options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry()){
      //Matrix_d tm(1,1);
      //AncestryProbs.resize( numCompositeLoci,tm );
      AncestryProbs = new Matrix_d[numCompositeLoci];
    }

    X_posn = 9999;
    string s1("\"X\"");
  // set size of locus ancestry array
    for( unsigned int j = 0; j < numChromosomes; j++ ){
        if( chrm[j]->GetLabel(0) != s1 ){// if not X chromosome, set number of elements to 2, num loci
            LocusAncestry(j).SetNumberOfElements( 2, chrm[j]->GetSize() );
            gametes.push_back(2);
        }
        else if( sex != 2 ){
            LocusAncestry(j).SetNumberOfElements( 1, chrm[j]->GetSize() );
            gametes.push_back(1);
            X_posn = j;
        }
        else{
            LocusAncestry(j).SetNumberOfElements( 2, chrm[j]->GetSize() );
            gametes.push_back(2);
            X_posn = j;
        }
        if( options->getPopulations() == 1 ) LocusAncestry(j).SetElements(0);

    }
    //retrieve genotypes
    vector<unsigned short> empty( 0, 0 );
    genotype.resize( numCompositeLoci, empty ); // stl vector of genotypes
    Data->GetGenotype(mynumber,options,Loci,&genotype);

    // loop over composite loci to set possible haplotype pairs compatible with genotype 
    for(int j=0;j<numCompositeLoci;++j) Loci(j)->SetPossibleHaplotypes(&PossibleHaplotypes[j],genotype[j]);
}

void Individual::InitialiseAffectedsOnlyScores(int L, int K){
  AffectedsScore.SetNumberOfElements(L, K);
  AffectedsVarScore.SetNumberOfElements(L, K);
  AffectedsInfo.SetNumberOfElements(L, K);
}
void Individual::InitialiseAncestryScores(int L, int K){
  AncestryScore.SetNumberOfElementsWithDimensions(L, 2 * K, 1);
  AncestryInfo.SetNumberOfElementsWithDimensions(L, 2 * K, 2 * K);
  AncestryVarScore.SetNumberOfElements(L, K);
  AncestryInfoCorrection.SetNumberOfElements(L, K);
  B.SetNumberOfElements(K, K);
}

Individual::~Individual()
{
  delete[] PossibleHaplotypes;
  //if(AncestryProbs)
    //TODO: delete this properly
    delete[] AncestryProbs;
   delete []f[0];
   delete []f[1];

}

vector< unsigned short >& Individual::getGenotype(unsigned int locus)
{
  return genotype[locus];
}

Vector_i Individual::getPossibleHaplotypes(unsigned int locus){
  return PossibleHaplotypes[locus];
}

//Indicates whether a genotype is missing at a locus
bool Individual::IsMissing(unsigned int locus)
{
  unsigned int count = 0;
  int NumberOfLoci = genotype[locus].size()/2;
  for(int i=0;i<NumberOfLoci;i++){
    count += genotype[locus][i*2]; // even-numbered elements
    //count += genotype[locus][i*2+1]; // odd-numbered elements
  }
  return (count == 0);
}

Matrix_d& Individual::getAdmixtureProps()
{
  return AdmixtureProps;
}

void Individual::setAdmixtureProps(Matrix_d a)
{
  AdmixtureProps = a;
}

Matrix_d& Individual::getAdmixturePropsX()
{
  return XAdmixtureProps;
}

void Individual::setAdmixturePropsX(Matrix_d a)
{
  XAdmixtureProps = a;
}

int Individual::getSex()
{
   return sex;
}

vector<bool>& Individual::getXi(unsigned int locus)
{
  return _xi[locus];
}

const vector< vector<bool> >&
Individual::getXi()
{
   return _xi;
}

Vector_i Individual::getSumXi()
{
   return sumxi;
}

double Individual::getSumrho0()
{
   return Sumrho0;
}

double Individual::getSumrho()
{
   double sumrho = 0;
   for( unsigned int g = 0; g < _rho.size(); g++ )
      sumrho += _rho[g];
   return sumrho;
}

vector<double> Individual::getRho()
{
   return _rho;
}

Vector_i Individual::GetLocusAncestry( int chrm, int locus )
{
   return LocusAncestry(chrm).GetColumn( locus );
}

//temporary function until ancestry ScoreTests are moved into here
Matrix_d Individual::getAncestryProbs(int locus)
{
  //this function used only for affectedonly and ancestry scoretests
    return( AncestryProbs[ locus ] );
}

double Individual::getLogPosteriorProb()
{
   return LogPosterior;
}

void Individual::UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates, Vector_d &poptheta, bool ModelIndicator,
Matrix_d *Covariates0)
{
  Vector_d avgtheta;
  if(ModelIndicator )
    avgtheta = AdmixtureProps.RowMean();
  else
    avgtheta = AdmixtureProps.GetColumn(0);
  for( int k = 0; k < Populations - 1; k++ )
    (*Covariates0)( i, NoCovariates - Populations + k + 1 )
      = avgtheta( k + 1 ) - poptheta( k + 1 );
}

// Metropolis update for admixture proportions theta, taking acceptance probability p as argument
void Individual::Accept_Reject_Theta( double p, bool xdata, int Populations, bool ModelIndicator )
{
  bool test = true;
  // loop over populations: if element of Dirichlet parameter vector is 0, do not update corresponding element of 
  // admixture proportion vector
  for( int k = 0; k < Populations; k++ ){
    if( (Theta)( k, 0 ) == 0.0 )
      test = false;
    else if( ModelIndicator && (Theta)( k, 1 ) == 0.0 )
      test = false;
  }

  // generic Metropolis rejection step
  if( p < 0 ){
     if( log(myrand()) < p && test ){
        setAdmixtureProps(Theta);
        if( xdata )
           setAdmixturePropsX(ThetaX);
     }
  }
  else{
     setAdmixtureProps(Theta);
     if( xdata )
        setAdmixturePropsX(ThetaX);
  }
}

// these two functions return ratio of likelihoods of new and old values of population admixture
// in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in. 
 
// should have just one function to get the likelihood in the regression model, given a value of population admixture
// should be generic method for GLM, given Xbeta, Y and probability distribution 
// then should calculate ratio in Metropolis step 
//  
double Individual::AcceptanceProbForTheta_LogReg( int i, int TI, bool ModelIndicator,int Populations,
					      int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, MatrixArray_d &ExpectedY, 
						  MatrixArray_d &Target, Vector_d &poptheta) 
{
  double prob, Xbeta = 0;
  // TI = Target Indicator, indicates which outcome var is used
  Vector_d avgtheta;
  // calculate mean of parental admixture proportions
  if( ModelIndicator )
    avgtheta = Theta.RowMean() - poptheta;
  else
    avgtheta = Theta.GetColumn(0) - poptheta;

  for( int jj = 0; jj < NoCovariates - Populations + 1; jj++ )
    Xbeta += Covariates0( i, jj ) * beta( TI )( jj ,0 );
  for( int k = 0; k < Populations - 1; k++ ){
    //? Old code had 0 instead of TI in for index of beta.
    Xbeta += avgtheta( k ) * beta(TI)( NoCovariates - Populations + k + 1, 0 );}
  double newExpectedY = 1 / ( 1. + exp( -Xbeta ) );
  if( Target( TI )( i, 0 ) == 1 )
    prob = newExpectedY / ExpectedY( TI )( i, 0 );
  else
    prob = ( 1 - newExpectedY ) / ( 1 - ExpectedY( TI )( i, 0 ) );

  return( log(prob) );
} 

double Individual::AcceptanceProbForTheta_LinearReg( int i, int TI,  bool ModelIndicator,int Populations,
						 int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, MatrixArray_d &ExpectedY,
						 MatrixArray_d &Target, Vector_d &poptheta, Vector_d &lambda)
{
  double prob, Xbeta = 0;
  Vector_d avgtheta;
  if( ModelIndicator )
    avgtheta = Theta.RowMean() - poptheta;
  else
    avgtheta = Theta.GetColumn(0) - poptheta;

  for( int jj = 0; jj < NoCovariates - Populations + 1; jj++ )
    Xbeta += Covariates0( i, jj ) * beta( TI )( jj, 0 );
  for( int k = 0; k < Populations - 1; k++ ){
    Xbeta += avgtheta( k ) * beta( TI )( NoCovariates - Populations + k + 1, 0 );
  }

  prob = 0.5 * lambda( TI ) * (( ExpectedY( TI )( i, 0 ) - Target( TI )( i, 0 ) ) * ( ExpectedY( TI )( i, 0 ) - Target( TI )( i, 0 ) )
			       - ( Xbeta - Target( TI )( i, 0 ) ) * ( Xbeta - Target( TI )( i, 0 ) ) );

  return( prob );
}

double Individual::AcceptanceProbForTheta_XChrm(std::vector<double> &sigma, int Populations )
{
   int gametes = 1;
   if( sex == 2 )
      gametes = 2;
   double p = 0;
   Matrix_d ThetaOld = AdmixtureProps, ThetaXOld = XAdmixtureProps;
   for( int g = 0; g < gametes; g++ ){
       p += gsl_sf_lngamma( sigma[g]*Theta.GetColumn(g).Sum() )
         - gsl_sf_lngamma( sigma[g]*ThetaOld.GetColumn(g).Sum() );
      for( int k = 0; k < Populations; k++ ){
         p += gsl_sf_lngamma( sigma[g]*ThetaOld(k,g) ) - gsl_sf_lngamma( sigma[g]*Theta(k,g) );
         p += (sigma[g]*Theta(k,g)-1.0)*log(ThetaX(k,g)) - (sigma[g]*ThetaOld(k,g)-1.0)*log(ThetaXOld(k,g));
      }
   }
   return p;
}

void Individual::SampleIndividualParameters( int i, Vector_d *SumLogTheta, AlleleFreqs *A, int iteration , MatrixArray_d *Target,
					     Vector_i &OutcomeType, MatrixArray_d &ExpectedY, Vector_d &lambda, int NoCovariates, 
					     Matrix_d &Covariates0,MatrixArray_d &beta, Vector_d &poptheta, 
					     AdmixOptions* options, Chromosome **chrm, 
					     vector<Vector_d> alpha, bool _symmetric, vector<bool> _admixed, double rhoalpha, 
					     double rhobeta,vector<double> sigma, double outcome, double DInvLink)
{
  double u;
  MyOutcome = outcome;

  if( options->getAnalysisTypeIndicator() > 1 ){
   // sample missing values of outcome variable
    for( int k = 0; k < Target->GetNumberOfElements(); k++ ){
      if( (*Target)(k).IsMissingValue( i, 0 ) ){
	if( !OutcomeType(k) ) // linear regression
	  (*Target)(k)( i, 0 ) = gennor( ExpectedY(k)( i, 0 ), 1 / sqrt( lambda(k) ) );
	else{// logistic regression
	  u = myrand();
	  if( u * ExpectedY(k)( i, 0 ) < 1 )
	    (*Target)(k)( i, 0 ) = 1;
	  else
	    (*Target)(k)( i, 0 ) = 0;
	}
      }
    }
  }

  //resets
   sumxi.SetElements( 0 );
 
   // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1  
    SumLocusAncestry.SetNumberOfElements(options->getPopulations(),2);
    if( options->getModelIndicator() ){//random mating model
      Theta.SetNumberOfElements( options->getPopulations(), 2 );
    }
    else{
      Theta.SetNumberOfElements( options->getPopulations(), 1 );
    }
    
    if(A->getLoci()->isX_data() ){
      if( sex == 1 ){
	ThetaX.SetNumberOfElements( options->getPopulations(), 1 );
	SumLocusAncestry_X.SetNumberOfElements( options->getPopulations(), 1 );
      }
      else{
	ThetaX.SetNumberOfElements( options->getPopulations(), 2 );
	SumLocusAncestry_X.SetNumberOfElements( options->getPopulations(), 2 );
      }
    }
    // next step should be to calculate transition matrices 

    // sampling locus ancestry requires calculation of forward probability vectors alpha in HMM 
    SampleLocusAncestry(chrm, options, A);
    
    // sample individual admixture proportions theta
    // should be modified to allow a population mixture component model   
    SampleParameters(i, options, A, chrm, alpha, _symmetric, _admixed,rhoalpha, rhobeta, iteration, sigma);

    double p = 0;
    if( options->getAnalysisTypeIndicator() == 2 && !options->getScoreTestIndicator() ){
      p = AcceptanceProbForTheta_LinearReg( i, 0, options->getModelIndicator(),options->getPopulations(),
					    NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta,lambda); 
    }
    else if( (options->getAnalysisTypeIndicator() == 3 || options->getAnalysisTypeIndicator() == 4) && !options->getScoreTestIndicator() ){
      p = AcceptanceProbForTheta_LogReg( i, 0, options->getModelIndicator(),options->getPopulations(),
					 NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta); 
    }
    else if( options->getAnalysisTypeIndicator() == 5 ){
      for( int k = 0; k < Target->GetNumberOfElements(); k++ ){
	if( OutcomeType( k ) )
	  p += AcceptanceProbForTheta_LogReg( i, k, options->getModelIndicator(), options->getPopulations(),
					      NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta); 
	else
	   p += AcceptanceProbForTheta_LinearReg( i, k, options->getModelIndicator(), options->getPopulations(),
						  NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta,lambda);
      }
    }
    if( A->getLoci()->isX_data() && !options->getXOnlyAnalysis() )
      p += AcceptanceProbForTheta_XChrm( sigma, options->getPopulations());
    
    Accept_Reject_Theta(p, A->getLoci()->isX_data(),options->getPopulations(), options->getModelIndicator() );

    if( options->getAnalysisTypeIndicator() > 1 )
      UpdateAdmixtureForRegression(i,options->getPopulations(), NoCovariates, poptheta, options->getModelIndicator(),&(Covariates0));
    for( int k = 0; k < options->getPopulations(); k++ ){
      (*SumLogTheta)( k ) += log( AdmixtureProps( k, 0 ) );
      if(options->getModelIndicator() && !options->getXOnlyAnalysis() )
	(*SumLogTheta)( k ) += log( AdmixtureProps( k, 1 ) );
    }
    
    // compute conditional distribution of locus ancestry if required for score tests
    if( options->getTestForAffectedsOnly()|| options->getTestForLinkageWithAncestry() )
      for( unsigned int j = 0; j < numChromosomes; j++ )
	for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ )
	  chrm[j]->getAncestryProbs( jj, &AncestryProbs[chrm[j]->GetLocus(jj)] );
    
    //Update affecteds only scores
    if(iteration > options->getBurnIn()){
      if( options->getAnalysisTypeIndicator() == 0 ) 
       UpdateScoreForLinkageAffectedsOnly(options->getPopulations(), options->getModelIndicator(), 
					  A->getLoci()->GetNumberOfCompositeLoci() );
      else if( options->getTestForAffectedsOnly() && MyOutcome == 1 ){
	UpdateScoreForLinkageAffectedsOnly(options->getPopulations(), options->getModelIndicator(), 
					   A->getLoci()->GetNumberOfCompositeLoci() );
      }
    }
    if(iteration > options->getBurnIn() && options->getTestForLinkageWithAncestry() ){
      Matrix_d Xcov(options->getPopulations(), 1);
      
      Xcov(options->getPopulations()-1, 0) = 1;
      for( int k = 0; k < options->getPopulations() - 1; k++ ){
	Xcov(k,0) = getAdmixtureProps()( k, 0 ); 
      }
      B += Xcov * Xcov.Transpose() * DInvLink;
    }
}

// Samples individual admixture proportions conditional on sampled values of ancestry at loci where 
// jump indicator xi is 1, population admixture distribution parameters alpha, and likelihood from regression 
// model (if there is one)  
// should have an alternative function to sample population mixture component membership and individual admixture proportions
// conditional on genotype, not sampled locus ancestry
void Individual::SampleParameters( int ind, AdmixOptions* options, AlleleFreqs *A, Chromosome **chrm, vector<Vector_d> alpha, 
			      bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta, int iteration, 
			      vector<double> sigma)
{


    double L = A->getLoci()->GetLengthOfGenome(), L_X=0.0;
    if( A->getLoci()->isX_data() ) L_X = A->getLoci()->GetLengthOfXchrm();
    
    //SumN is the number of arrivals between each pair of adjacent loci  
    unsigned int SumN[] = {0,0};
    unsigned int SumN_X[] = {0,0};
    
    // sample sum of intensities parameter rho
    if( options->getRhoIndicator() ){
      
      SampleNumberOfArrivals(A, options, chrm, SumN, SumN_X); // SumN is number of arrivals
      
      SampleRho( options->getXOnlyAnalysis(), options->getModelIndicator(), A->getLoci()->isX_data(), rhoalpha, rhobeta, L, L_X, 
		 SumN, SumN_X);
    }
    
    Vector_d vectemp;//used to hold sample from theta posterior
    // if no regression model, sample admixture proportions theta as a conjugate Dirichlet posterior   
    if( options->getXOnlyAnalysis() ){
      vectemp = gendirichlet( alpha[0] + SumLocusAncestry_X.GetColumn(0) );
      Theta.SetColumn( 0, vectemp );
    }
    else if( options->getModelIndicator() ){//random mating model
      for( unsigned int g = 0; g < 2; g++ ){
	if( options->getAnalysisTypeIndicator() > -1 )
	  vectemp = gendirichlet( alpha[0] + SumLocusAncestry.GetColumn(g) );
	else
	  vectemp = gendirichlet( alpha[g] + SumLocusAncestry.GetColumn(g) );
	Theta.SetColumn( g, vectemp );
      }
      if( A->getLoci()->isX_data() ){//nonrandom mating model
	for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	  vectemp = gendirichlet( Theta.GetColumn(g)*sigma[g]
				  + SumLocusAncestry_X.GetColumn(g) );
	  ThetaX.SetColumn( g, vectemp );
	}
      }
    }
    else{// no random mating model
      vectemp = gendirichlet( alpha[0] + SumLocusAncestry.RowSum() );
      Theta.SetColumn( 0, vectemp );
    }

    // calculate posterior density conditional on realized locus ancestry states, jump indicators and num arrivals 
    if( options->getMLIndicator() && ind == 0 && iteration > options->getBurnIn() )
      CalculateLogPosterior(options,A->getLoci()->isX_data(), alpha, _symmetric,
			    _admixed,rhoalpha, rhobeta, L, L_X, SumN, SumN_X);
    
}

void Individual::SampleLocusAncestry(Chromosome **chrm, AdmixOptions *options, AlleleFreqs *A){
  // Loops over loci to sample locus ancestry and jump indicators xi 
  //updates allele counts in A
  int locus;
  double q, Prob;

  //fill f with values in AlleleFreqs;
  A->getLociCorrSummary(f);
  
  Sumrho0 = 0;
  locus = 0;
  // f0 and f1 are arrays of scalars of the form exp - rho*x, where x is distance between loci
  // required to calculate transition matrices 
  if( options->getRhoIndicator() ){
    for( unsigned int j = 0; j < numChromosomes; j++ ){
      locus++;
      for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
	f[0][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho[0] );
	if( options->getModelIndicator() ){
	  f[1][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho[1] );
	}
	else
	  f[1][locus] = f[0][locus];
	locus++;
      }
    }
  }
 
  locus = 0;
  bool isdiploid;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    locus = chrm[j]->GetLocus(0);

    //Update Forward/Backward probs in HMM
    if( j != X_posn ){
      chrm[j]->UpdateParameters( this, A,AdmixtureProps, options, f, false, true );
      isdiploid = true;
   }
    else if( options->getXOnlyAnalysis() ){
      chrm[j]->UpdateParameters( this, A,AdmixtureProps, options, f, false, false );
      isdiploid = false;
    }
    else if( sex == 1 ){
      chrm[j]->UpdateParameters( this, A,XAdmixtureProps, options, f, false, false );
      isdiploid = false;
    }
    else{
      chrm[j]->UpdateParameters( this, A,XAdmixtureProps, options, f, false, true );
      isdiploid = true;
    }

    //Sample locus ancestry
    chrm[j]->SampleForLocusAncestry(&LocusAncestry(j),isdiploid);

    //loop over loci on current chromosome and update allele counts and Ancestry probs if needed
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      int loc =  chrm[j]->GetLocus(jj);
      if( !(IsMissing(j)) ){
	  if(isdiploid)
	    A->UpdateAlleleCounts( loc, PossibleHaplotypes[loc], LocusAncestry(j).GetColumn(jj) );
	  else
	    A->UpdateAlleleCounts_HaploidData( loc, genotype[loc], LocusAncestry(j)(0,jj) );
	}
     }   

    //locus++;
    //update jump indicators xi  
    for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
      locus = chrm[j]->GetLocus(jj);    
      for( unsigned int g = 0; g < gametes[j]; g++ ){
	if( LocusAncestry(j)(g,jj-1) == LocusAncestry(j)(g,jj) ){
	  if( options->getModelIndicator() && g == 1 ){
	    q = AdmixtureProps( LocusAncestry(j)(g,jj), g );
	  } else {
	    q = AdmixtureProps( LocusAncestry(j)(g,jj), 0 );
	  }
	  Prob = (1 - f[g][locus]) / ( 1 - f[g][locus] + f[g][locus] / q );
	  if( Prob > myrand() ){
	    _xi[g][locus] = true;
	    sumxi(locus)++;
	  } else {
	    _xi[g][locus] = false;
	    Sumrho0 += A->getLoci()->GetDistance( locus );
	  }
	} else {
	  _xi[g][locus] = true;
	  sumxi(locus)++;
	}
      }
      //locus++;
    }
    locus = chrm[j]->GetLocus(0);
    // sum ancestry states over loci where jump indicator is 1
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      for( unsigned int g = 0; g < gametes[j]; g++ ){
	if( _xi[g][locus] ){
	  if( j != X_posn )
	    SumLocusAncestry( LocusAncestry(j)( g, jj ), g )++;
	  else
	    SumLocusAncestry_X( LocusAncestry(j)( g, jj ), g )++;
	}
      }
      locus++;
    }
  }//end chromosome loop

}

void Individual::SampleNumberOfArrivals(AlleleFreqs *A, AdmixOptions *options, Chromosome **chrm, 
					unsigned int SumN[], unsigned int SumN_X[]){
  // samples number SumN of arrivals between each pair of adjacent loci, 
  // conditional on jump indicators xi and sum of intensities rho
  // total number SumN is used for conjugate update of sum of intensities 
  double q, rho = 0.0;
  int locus = 0;
  int ran = 0;
  if( myrand() < 0.5 ) ran = 1;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    locus++;
    for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
      double delta = A->getLoci()->GetDistance(locus);
      for( unsigned int g = 0; g < gametes[j]; g++ ){
	if( _xi[g][locus] ){
	  if( options->getXOnlyAnalysis() ){
	    rho = _rho[0];
	    q = AdmixtureProps( LocusAncestry(j)(0,jj-ran), 0 );
	  }
	  else if( options->getModelIndicator() && g == 1 ){
	    q = AdmixtureProps( LocusAncestry(j)(g,jj-ran), g );
                    if( j != X_posn )
		      rho = _rho[g];
                    else
		      rho = _rho_X[g];
	  } else {
	    q = AdmixtureProps( LocusAncestry(j)(g,jj-ran), 0 );
	    if( j != X_posn )
	      rho = _rho[0];
	    else
	      rho = _rho_X[0];
	  }
	  double u = myrand();
	  //                 double deltadash = -log( 1 - u*( 1 - exp(-q*rho*delta) ) ) / (q*rho);
	  double deltadash = -log( 1 - u*( 1 - exp(-rho*delta) ) ) / (rho);
	  unsigned int sample = genpoi( rho*(delta - deltadash) );
	  if( j != X_posn )
	    SumN[g] += sample + 1;
	  else
	    SumN_X[g] += sample + 1;
	}
      }
      locus++;
    }
  }
}

void Individual::SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta, double L, double L_X, 
			   unsigned int SumN[], unsigned int SumN_X[]){
  // Samples sum of intensities parameter as conjugate gamma with Poisson likelihood
  // SumN is the number of arrivals between each pair of adjacent loci
  if(XOnly ){
    do{
      _rho[0] = gengam( rhobeta + L_X, rhoalpha + (double)SumN_X[0] );
    }while( _rho[0] > TruncationPt || _rho[0] < 1.0 );
  }
  else if(RandomMatingModel ){
    for( unsigned int g = 0; g < 2; g++ ){
      do{
	_rho[g] = gengam( rhobeta + L, rhoalpha + (double)SumN[g] );
      }while( _rho[g] > TruncationPt || _rho[g] < 1.0 );
    }
    if(X_data  ){
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	do{
	  _rho_X[g] = gengam( rhobeta + L_X, rhoalpha + (double)SumN_X[g] );
	}while( _rho[g] > TruncationPt || _rho[g] < 1.0 );
      }
    }
  }
  else{
    _rho[0] = gengam( rhobeta + 2*L, rhoalpha + (double)(SumN[0] + SumN[1]) );
  }
}

void Individual::ResetScores(AdmixOptions *options){
  if( options->getTestForAffectedsOnly() ){
    AffectedsScore.SetElements(0);
    AffectedsVarScore.SetElements(0);
    AffectedsInfo.SetElements(0);
  }
  if( options->getTestForLinkageWithAncestry() ){
    AncestryScore.SetElements(0);
    AncestryInfo.SetElements(0);
    AncestryInfoCorrection.SetElements(0);
    AncestryVarScore.SetElements(0);
    B.SetElements(0);    
  }
}

void Individual::UpdateScoreForLinkageAffectedsOnly(int Populations, bool ModelIndicator, int NumberOfCompositeLoci){
  // Different from the notation in McKeigue et  al. (2000). McKeigue
  // uses P0, P1, P2, which relate to individual admixture as follows;
  // P0 = ( 1 - theta0 ) * ( 1 - theta1 )
  // P1 = ( 1 - theta0 ) * theta1 + theta0 * ( 1 - theta1 )
  // P2 = theta0 * theta1

  // Test in McKeigue et al. (2000) was on r defined on the positive
  // real line.  This test is on log r, which should have better
  // asymptotic properties.

  double theta[2];//paternal and maternal admixture proportions
  //Matrix_d AncestryProbs;//conditional locus ancestry probs

  for( int k = 0; k < Populations; k++ ){
    theta[0] = AdmixtureProps( k, 0 );
    if( ModelIndicator )
      theta[1] = AdmixtureProps( k, 1 );
    else
      theta[1] = AdmixtureProps( k, 0 );

    for(int locus = 0; locus< NumberOfCompositeLoci; ++locus){
      //AncestryProbs = getAncestryProbs(locus);
      AffectedsScore(locus,k)+= 0.5*( AncestryProbs[locus](k,1) + 2.0*AncestryProbs[locus](k,2) - theta[0] - theta[1] );
      AffectedsVarScore(locus, k)+= 0.25 *( AncestryProbs[locus](k,1)*(1.0 - AncestryProbs[locus](k,1)) + 
					    4.0*AncestryProbs[locus](k,2)*AncestryProbs[locus](k,0)); 
      AffectedsInfo(locus, k)+= 0.25* ( theta[0]*( 1.0 - theta[0] ) + theta[1]*( 1.0 - theta[1] ) );
    }
  }
}

void Individual::UpdateScoreForAncestry(double phi, double YMinusEY, double DInvLink, int NumberOfCompositeLoci,int Populations)
{
  //Updates score stats for test for association with locus ancestry
  //now use Rao-Blackwellized estimator by replacing realized ancestries with their expectations
  //Notes: 1/phi is dispersion parameter
  //       = lambda(0) for linear regression, = 1 for logistic
  //       YMinusEY = Y - E(Y) = Y - g^{-1}(\eta_i)
  //       VarX = Var(X)
  //       DInvLink = {d  g^{-1}(\eta)} / d\eta = derivative of inverse-link function
  //Xcov is a vector of covariates
  //Note that only the intercept and admixture proportions are used.
  // X is (A, Xcov)'  

  int dim = 2 * Populations;
  Matrix_d Aprobs;
  double xBx;
  Matrix_d X, Xcov;
  Vector_d VarA, temp;
  //double YMinusEY = MyOutcome - EY;

 
  X.SetNumberOfElements(dim, 1);
  Xcov.SetNumberOfElements(Populations, 1);
  VarA.SetNumberOfElements(Populations);
  temp.SetNumberOfElements(Populations);
      
      X( dim - 1, 0 ) = 1;//intercept
      Xcov(Populations-1, 0) = 1;
      //set covariates 
      for( int k = 0; k < Populations - 1; k++ ){
	X( k + Populations, 0 ) = getAdmixtureProps()( k, 0 );
	Xcov(k,0) = getAdmixtureProps()( k, 0 ); 
      }
      
      for(int locus =0; locus < NumberOfCompositeLoci; ++locus){
	
	Aprobs = getAncestryProbs(locus);//conditional locus ancestry probs      
	
	for( int k = 0; k < Populations ; k++ ){
	  X(k,0) = Aprobs(k,1) + 2.0 * Aprobs(k,2);//Conditional expectation of ancestry
	  VarA(k) = Aprobs(k,1)*(1.0 - Aprobs(k,1)) + 4.0*Aprobs(k,2)*Aprobs(k,0);//conditional variances
	}
	
	temp =  Xcov.GetColumn(0);
	HH_svx(B, &temp);
	xBx = (Xcov.Transpose() * temp.ColumnMatrix())(0,0);
	
	AncestryScore(locus) += X * YMinusEY * phi;
	AncestryInfo(locus) += (X * X.Transpose()) * DInvLink * phi;
	
	for( int k = 0; k < Populations ; k++ ){
	  AncestryInfoCorrection(locus,k) += VarA(k) * (DInvLink *phi - phi * phi * DInvLink * DInvLink * xBx) ;   
	  AncestryVarScore(locus,k) += VarA(k) * phi * phi * YMinusEY * YMinusEY;
	}
      }//end locus loop
}

void Individual::SumScoresForLinkageAffectedsOnly(int j,int Populations, Matrix_d *SumAffectedsScore, 
				      Matrix_d *SumAffectedsVarScore,Matrix_d *SumAffectedsScore2, Matrix_d *SumAffectedsInfo){
  for( int kk = 0; kk < Populations; kk++ ){
    (*SumAffectedsScore)(j,kk) += AffectedsScore(j,kk);
    (*SumAffectedsVarScore)(j,kk) += AffectedsVarScore(j,kk);
    (*SumAffectedsInfo)(j,kk) += AffectedsInfo(j,kk);
    (*SumAffectedsScore2)(j,kk) +=  AffectedsScore(j, kk) * AffectedsScore(j, kk);
  }
}

void Individual::SumScoresForAncestry(int j, int Populations, Matrix_d *score, Matrix_d *info, 
				      Matrix_d *SumAncestryScore, Matrix_d *SumAncestryInfo, Matrix_d *SumAncestryScore2,
				      Matrix_d *SumAncestryVarScore){
  CentredGaussianConditional(Populations,AncestryScore(j), AncestryInfo(j), score, info );
  //accumulate over iterations     
  for( int k = 0; k < Populations ; k++ ){
    (*SumAncestryScore)(j,k) += (*score)(k,0);
    (*SumAncestryInfo)(j,k)  += (*info)(k,k) + AncestryInfoCorrection(j,k);
    (*SumAncestryScore2)(j,k) += (*score)(k,0) * (*score)(k,0);
    (*SumAncestryVarScore)(j,k) += AncestryVarScore(j,k);
  }
}

// unnecessary duplication of code - should use same method as for > 1 population
void Individual::OnePopulationUpdate( int i, MatrixArray_d *Target, Vector_i &OutcomeType, MatrixArray_d &ExpectedY, Vector_d &lambda, 
				     int AnalysisTypeIndicator )
{
  Vector_i ancestry(2);
  for( int k = 0; k < Target->GetNumberOfElements(); k++ ){
    if( AnalysisTypeIndicator > 1 ){
      if( (*Target)(k).IsMissingValue( i, 0 ) ){
	if( !OutcomeType(k) )
	  (*Target)(k)( i, 0 ) = gennor( ExpectedY(k)( i, 0 ), 1 / sqrt( lambda(k) ) );
	else{
	  if( myrand() * ExpectedY(k)( i, 0 ) < 1 )
	    (*Target)(k)( i, 0 ) = 1;
	  else
	    (*Target)(k)( i, 0 ) = 0;
	}
      }
    }
  }
  // sampled alleles should be stored in Individual objects, then summed over individuals to get counts
  // is this extra update necessary? isn't this method called anyway, irrespective of whether there is only one population?        
 //  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
//     A->UpdateAlleleCounts(j,getPossibleHaplotypes(j), ancestry );
//   }
}

void Individual::InitializeChib(Matrix_d theta, Matrix_d thetaX, vector<double> rho, vector<double> rhoX, 
				AdmixOptions *options, AlleleFreqs *A, Chromosome **chrm, double rhoalpha, double rhobeta, 
				vector<Vector_d> alpha, vector<bool> _admixed, chib *MargLikelihood, std::ofstream *LogFileStreamPtr)
//Computes LogPrior and LogLikelihood used for Chib Algorithm
{
   double LogPrior=0, LogLikelihoodAtEst;
   *LogFileStreamPtr << "Calculating posterior at individual admixture\n"
                    << theta << "and rho\n" << rho[0] << " " << rho[1] << endl;
   if( options->getXOnlyAnalysis() ){
     LogLikelihoodAtEst = getLogLikelihoodXOnly( options, A, chrm, theta, rho );
      if( options->getRho() == 99 ){
         LogPrior = -log( options->getTruncPt() - 1.0 );
      }
      else if( options->getRho() == 98 ){
         LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) );
      }
      else{
         LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] );
         LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
      }
      LogPrior += getDirichletLogDensity( alpha[0], theta.GetColumn(0) );
   }
   else if( A->getLoci()->isX_data() ){
     LogLikelihoodAtEst = getLogLikelihood( options, A, chrm, theta, rho, thetaX, rhoX );
      if( options->getRho() == 99 ){
         LogPrior = -4.0*log( options->getTruncPt() - 1.0 );
      }
      else if( options->getRho() == 98 ){
         LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) )
            -log( rho[1]*(log( options->getTruncPt() ) ) )
            -log( rhoX[0]*(log( options->getTruncPt() ) ) )
            -log( rhoX[1]*(log( options->getTruncPt() ) ) );
      }
      else{
         LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] )
            + getGammaLogDensity( rhoalpha, rhobeta, rhoX[0] )
            + getGammaLogDensity( rhoalpha, rhobeta, rho[1] )
            + getGammaLogDensity( rhoalpha, rhobeta, rhoX[1] );
         LogPrior /= gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0);
      }
      LogPrior += getDirichletLogDensity( alpha[0], theta.GetColumn(0) )
         + getDirichletLogDensity( alpha[0], thetaX.GetColumn(0) )
         + getDirichletLogDensity( alpha[1], theta.GetColumn(1) )
         + getDirichletLogDensity( alpha[1], thetaX.GetColumn(1) );
      LogLikelihoodAtEst = getLogLikelihood( options, A, chrm, theta, rho, thetaX, rhoX );
   }
   else{
      if( options->getPopulations() > 1 ){
         LogLikelihoodAtEst = getLogLikelihood( options, A, chrm, theta, rho, thetaX, rhoX );
         if( _admixed[0] ){
            if( options->getRho() == 99 ){
               LogPrior = -log( options->getTruncPt() - 1.0 );
            }
            else if( options->getRho() == 98 ){
               LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) );
            }
            else{
               LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] );
               LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
            }
            LogPrior += getDirichletLogDensity( alpha[0], theta.GetColumn(0) );
         }
         if( _admixed[1] ){

            if( options->getRho() == 99 ){
               LogPrior -= log( options->getTruncPt() - 1.0 );
            }
            else if( options->getRho() == 98 ){
               LogPrior -= log( rho[1]*(log( options->getTruncPt() ) ) );
            }
            else{
               LogPrior += getGammaLogDensity( rhoalpha, rhobeta, rho[1] );
               LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
            }
            LogPrior += getDirichletLogDensity( alpha[1], theta.GetColumn(1) );
         }
      }
      else{
         LogLikelihoodAtEst = getLogLikelihoodOnePop(A);
      }
   }
   if( A->IsRandom() ){
      for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
         for( int k = 0; k < options->getPopulations(); k++ ){
	   //CompositeLocus *locus = (CompositeLocus*)(*Loci)(j);
	   LogPrior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
         }
      }
   }
   MargLikelihood->setLogPrior( LogPrior );
   MargLikelihood->setLogLikelihood( LogLikelihoodAtEst );   
}

 // this function computes marginal likelihood by the Chib algorithm.  Can be replaced with 
  // more efficient algorithm based on HMM likelihood
// Chib method for marginal likelihood should be rewritten to use the HMM likelihood, without sampling locus ancestry or arrivals
void Individual::ChibLikelihood(int i,int iteration, double *LogLikelihood, double *SumLogLikelihood, vector<double> MaxLogLikelihood, 
				AdmixOptions *options, Chromosome **chrm, vector<Vector_d> alpha, 
				vector<bool> _admixed, double rhoalpha, double rhobeta, MatrixArray_d &thetahat, 
				MatrixArray_d &thetahatX, vector<vector<double> > &rhohat, 
				vector<vector<double> > &rhohatX,std::ofstream *LogFileStreamPtr, chib *MargLikelihood, AlleleFreqs* A){
            
//           if( iteration <= options->getBurnIn() ){

  *LogLikelihood = getLogLikelihood(options, A, chrm);
    if( options->getPopulations() > 1 ){
      if( options->getRho() < 90 ){
	if( _admixed[0] ){
	  *LogLikelihood+=getGammaLogDensity( rhoalpha, rhobeta, _rho[0] );}
	if( _admixed[1] )
	  *LogLikelihood+=getGammaLogDensity( rhoalpha, rhobeta, _rho[1] );
      }
      else if( options->getRho() == 98 ){
	if( _admixed[0] )
	  *LogLikelihood -= log( _rho[0] );
	if( _admixed[1] )
	  *LogLikelihood -= log( _rho[1] );
      }
      *LogLikelihood+=
	getDirichletLogDensity(alpha[0],
			       getAdmixtureProps().GetColumn(0))
	+getDirichletLogDensity(alpha[1],
				getAdmixtureProps().GetColumn(1));
      if( *LogLikelihood > MaxLogLikelihood[i] ){
	*LogFileStreamPtr << getAdmixtureProps()
			 << _rho[0] << " " << _rho[1]
			 << endl << *LogLikelihood << endl
			 << iteration << endl;
	MaxLogLikelihood[i] = *LogLikelihood;
	if( iteration <= options->getBurnIn() ){
	  thetahat(i) = getAdmixtureProps();
	  rhohat[i] = _rho;
	  A->setAlleleFreqsMAP();
	  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
	    CompositeLocus *locus = (CompositeLocus*)(A->getLocus(j));
	    //locus->setAlleleFreqsMAP();
	    if( locus->GetNumberOfLoci() > 2 )
	      locus->setHaplotypeProbsMAP();
	  }
	  if( A->getLoci()->isX_data() ){
	    thetahatX(i) = getAdmixturePropsX();
	    rhohatX[i] = _rho_X;
	  }
	}
      }
    }
    else{//populations <=0
      if( *LogLikelihood > MaxLogLikelihood[i] ){

	MaxLogLikelihood[i] = *LogLikelihood;
	if( iteration <= options->getBurnIn() ){
	  A->setAlleleFreqsMAP();
	  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
	    CompositeLocus *locus = (CompositeLocus*)(A->getLocus(j));
	    //locus->setAlleleFreqsMAP();
	    if( locus->GetNumberOfLoci() > 2 )
	      locus->setHaplotypeProbsMAP();
	  }
	}
      }
    }
    if( options->getAnalysisTypeIndicator() == -1 ){
      if( iteration == options->getBurnIn() ){
	InitializeChib(thetahat(0), thetahatX(0), rhohat[0], rhohatX[0], 
		       options, A, chrm, rhoalpha, rhobeta, 
		       alpha, _admixed, MargLikelihood, LogFileStreamPtr);
      }
      if( iteration > options->getBurnIn() ){
	double LogPosterior = 0;
	if( options->getPopulations() > 1 )
	  LogPosterior = getLogPosteriorProb();
	if( A->IsRandom() ){
	  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
	    for( int k = 0; k < options->getPopulations(); k++ ){
	      //CompositeLocus *locus = (CompositeLocus*)(*Loci)(j);
	      LogPosterior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j,k) + A->GetAlleleCounts(j,k), 
						      A->getAlleleFreqsMAP(j, k) );
	    }
	  }
	}
	MargLikelihood->addLogPosteriorObs( LogPosterior );
	*SumLogLikelihood += *LogLikelihood;
      }
    }

}

double
Individual::getLogLikelihoodXOnly( AdmixOptions* options, AlleleFreqs *A, Chromosome **chrm, Matrix_d ancestry, vector<double> rho )
{
   double LogLikelihood = 0.0;
   _rhoHat = rho;
   AdmixtureHat = ancestry;
   Vector_d ff( A->GetNumberOfCompositeLoci() );

   for(int j=0;j<A->GetNumberOfCompositeLoci();++j)f[0][j] = ff(j);
  
   for( unsigned int jj = 1; jj < chrm[0]->GetSize(); jj++ ){
     f[0][jj] = exp( -A->getLoci()->GetDistance( jj ) * _rhoHat[0] );
   }
   chrm[0]->UpdateParameters( this, A,AdmixtureHat, options, f, true, false );
   LogLikelihood += chrm[0]->getLogLikelihood();

   return LogLikelihood;
}

//do we need 2 getLogLikelihood functions?
//This one is called by InitializeChib, called in turn by ChibLikelihood
//only used to compute marginal likelihood   
double
Individual::getLogLikelihood( AdmixOptions* options, AlleleFreqs *A, Chromosome **chrm, Matrix_d ancestry, vector<double> rho, Matrix_d ancestry_X, vector<double> rho_X )
{
  //Why are the forward, backward probs updated here when this has been done already?
   int locus = 0;
   _rhoHat = rho;
   AdmixtureHat = ancestry;
   double LogLikelihood = 0.0;
   Vector_d ff( A->GetNumberOfCompositeLoci() );

   for(int j=0;j<2;++j){
     for(int k;k<A->GetNumberOfCompositeLoci();++k)f[j][k] = ff(k);
   }
   
   for( unsigned int j = 0; j < numChromosomes; j++ ){      
      locus++;
      if( j != X_posn ){
	for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
	   f[0][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat[0] );
            if( options->getModelIndicator() ){
	      f[1][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat[1] );
            }
            else
               f[1][locus] = f[0][locus];
            locus++;
         }
         chrm[j]->UpdateParameters( this,A, AdmixtureHat, options, f, true, true);
         //LogLikelihood += chrm[j]->getLogLikelihood();
      }
      else{
         _rhoHat_X = rho_X;
         XAdmixtureHat = ancestry_X;
         for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
	   f[0][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat_X[0] );
            if( sex == 2 ){
	      f[1][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat_X[1] );
            }
            locus++;
         }
         if( sex == 1 )
	   chrm[j]->UpdateParameters( this,A, XAdmixtureHat, options, f, true, false);
         else
	   chrm[j]->UpdateParameters( this, A,XAdmixtureHat, options, f, true, true);
         //LogLikelihood += chrm[j]->getLogLikelihood();
      }
      LogLikelihood += chrm[j]->getLogLikelihood();
   }

   return LogLikelihood;
}

double
Individual::getLogLikelihoodOnePop(AlleleFreqs *A )
{
   double Likelihood = 0.0;
   Matrix_d Prob;
   for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
     if(!IsMissing(j)){
       A->GetGenotypeProbs(&Prob, j, genotype[j], getPossibleHaplotypes(j), true, true );
       Likelihood += log( Prob(0,0) );
     }
   }
   return Likelihood;
}

//do we need 2 getLogLikelihood functions?
//called at top of ChibLikelihood, used to compute marginal likelihood
double
Individual::getLogLikelihood( AdmixOptions* options, AlleleFreqs* A, Chromosome **chrm )
{
   int locus = 0;
   double LogLikelihood = 0.0;
   Vector_d ff( A->GetNumberOfCompositeLoci() );
 
   for(int j=0;j<2;++j){
     for(int k;k<A->GetNumberOfCompositeLoci();++k)f[j][k] = ff(k);
   }
   
   for( unsigned int j = 0; j < numChromosomes; j++ ){      
      locus++;
      if( j != X_posn ){
	for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
	  f[0][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho[0] );
            if( options->getModelIndicator() ){
	      f[1][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho[1] );
            }
            else
               f[1][locus] = f[0][locus];
            locus++;
         }
	chrm[j]->UpdateParameters( this, A,AdmixtureProps, options, f, false, true);
	//LogLikelihood += chrm[j]->getLogLikelihood();
      }
      else if( options->getXOnlyAnalysis() ){
	for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
	  f[0][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho[0] );
            locus++;
         }
         chrm[j]->UpdateParameters( this, A,AdmixtureProps, options, f, false, false);
         //LogLikelihood += chrm[j]->getLogLikelihood();
      }
      else{
	for( unsigned int jj = 1; jj < chrm[j]->GetSize(); jj++ ){
	   f[0][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho_X[0] );
            if( sex == 2 ){
	      f[1][locus] = exp( -A->getLoci()->GetDistance( locus ) * _rho_X[1] );
            }
            locus++;
         }
         if( sex == 1 )
	   chrm[j]->UpdateParameters( this, A,AdmixtureProps, options, f, false, false );
         else
	   chrm[j]->UpdateParameters( this, A,AdmixtureProps, options, f, false, true );
         //LogLikelihood += chrm[j]->getLogLikelihood();
      }
      LogLikelihood += chrm[j]->getLogLikelihood();
   }

   return LogLikelihood;
}

void Individual::CalculateLogPosterior(AdmixOptions *options, bool isX_data, vector<Vector_d> alpha, 
						 bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta, double L, 
				       double L_X, unsigned int SumN[],unsigned int SumN_X[]){

  LogPosterior = 0.0; 
  double IntConst1;
 {
    Vector_d alphaparams1, alphaparams0;
    if( options->getXOnlyAnalysis() ){
      LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[0],
					  rhobeta + L_X, _rhoHat[0] );
      if( options->getRho() > 90.0 )
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[0], rhobeta+L_X, 1.0, options->getTruncPt() );
      else
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[0], 1.0);
      LogPosterior -= log(IntConst1);
      alphaparams0 = alpha[0] + SumLocusAncestry_X.GetColumn(0);
      LogPosterior += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(0));
    }
    else if( isX_data ){
      for( unsigned int g = 0; g < 2; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[g],
					    rhobeta + L, _rhoHat[g] );
	if( options->getRho() > 90.0 )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[g], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[g], 1.0);
	LogPosterior -= log(IntConst1);
	alphaparams0 = alpha[g] + SumLocusAncestry.GetColumn(g);
	LogPosterior += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(g));
      }
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[g],
					    rhobeta + L_X, _rhoHat_X[g] );
	if( options->getRho() > 90.0 )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[g], rhobeta+L_X, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[g], 1.0);
	LogPosterior -= log(IntConst1);
	alphaparams0 = alpha[g] + SumLocusAncestry_X.GetColumn(g);
	LogPosterior += getDirichletLogDensity(alphaparams0,XAdmixtureHat.GetColumn(g));
      }
    }
    else if( _symmetric ){
      vector<double> x(2,0.0);
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[0],
				  rhobeta + L, _rhoHat[0] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[1],
				  rhobeta + L, _rhoHat[0] );
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[1],
				  rhobeta + L, _rhoHat[1] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[0],
				  rhobeta + L, _rhoHat[1] );
      double IntConst1, IntConst2;
      if( options->getRho() > 90.0 ){
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	IntConst2 = IntegratingConst(rhoalpha+(double)SumN[1], rhobeta+L, 1.0, options->getTruncPt() );
      }
      else{
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
	IntConst2 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
      }
      x[0] -= log(IntConst1);
      x[1] -= log(IntConst2);
      alphaparams0 = alpha[0] + SumLocusAncestry.GetColumn(0);
      x[0] += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(0));
      x[1] += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(1));
      alphaparams1 = alpha[1] + SumLocusAncestry.GetColumn(1);
      x[0] += getDirichletLogDensity(alphaparams1,AdmixtureHat.GetColumn(1));
      x[1] += getDirichletLogDensity(alphaparams1,AdmixtureHat.GetColumn(0));
      if( isnan(x[0]) || isinf(x[0]) )
	LogPosterior = x[1] - log(2.0);
      else if( isnan(x[1]) || isinf(x[1])  )
	LogPosterior = x[0] - log(2.0);
      else if( x[0] < x[1] )
	LogPosterior = x[1] + log( 1 + exp( x[0] - x[1] ) ) - log(2.0);
      else
	LogPosterior = x[0] + log( 1 + exp( x[1] - x[0] ) ) - log(2.0);
      if( isnan(LogPosterior) || isinf(LogPosterior) ){
	PR(alphaparams0);
	PR(alphaparams1);
	PR(AdmixtureHat.GetColumn(0));
           PR(AdmixtureHat.GetColumn(1));
           PR(x[0]);
           PR(x[1]);
           exit(0);
      }
    }
    else{
      if( _admixed[0] ){
	LogPosterior = getGammaLogDensity( rhoalpha + (double)SumN[0],
					   rhobeta + L, _rhoHat[0] );
	if( options->getRho() > 90.0 )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
           LogPosterior -= log( IntConst1 );
           alphaparams0 = alpha[0] + SumLocusAncestry.GetColumn(0);
           LogPosterior+=getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(0));
      }
      if( _admixed[1] ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[1],
					    rhobeta + L, _rhoHat[1] );
	if( options->getRho() > 90.0 )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[1], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
	LogPosterior -= log( IntConst1 );
	alphaparams1 = alpha[1] + SumLocusAncestry.GetColumn(1);
	LogPosterior+=getDirichletLogDensity(alphaparams1,AdmixtureHat.GetColumn(1));
      }
    }
  }
}

double Individual::IntegratingConst( double alpha, double beta, double a, double b )
{
   double I = gsl_cdf_gamma_P( b*beta, alpha, 1 ) - gsl_cdf_gamma_P( a*beta, alpha, 1);
   return I;
}
