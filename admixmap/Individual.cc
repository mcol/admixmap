#include "Individual.h"
#include "StringConvertor.h"

#define PR(x) cout << #x << " = " << x << endl;

Individual::Individual()
{
}

Individual::Individual(AdmixOptions* options, const Vector_s& data, Genome& Loci,Chromosome **chrm)
{
  //numChromosomes  = chrm.size();
  numChromosomes = Loci.GetNumberOfChromosomes();
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

    int TotalLoci = 0;
    //Can we not just count the length of Genome?
    //TotalLoci = Loci.GetNumberOfCompositeLoci();
    numCompLoci.resize( numChromosomes );
    for( unsigned int j = 0; j < numChromosomes; j++ ){
        numCompLoci[j] = chrm[j]->GetSize();
        for( unsigned int jj = 0; jj < numCompLoci[j]; jj++ ){
            int compLocus = chrm[j]->GetLocus(jj);
            TotalLoci += Loci(compLocus)->GetNumberOfLoci();
        }
    }

    if (options->IsPedFile() == 1) {
        if (data.size() != 2*TotalLoci + 1 + options->genotypesSexColumn()) {
            cout << "Error in formatting of line" << endl;
            exit(0);
        }
    } else {
        if (data.size() != TotalLoci + 1 + options->genotypesSexColumn()) {
            cout << "Error in formatting of line" << endl;
            exit(0);
        }
    }

    // Read sex value if present.
    if (options->genotypesSexColumn() == 1) {
        sex = StringConvertor::toInt(data[1]);
        if (sex > 2) {
            cout << "Error: sex must be coded as 0 - missing, 1 - male or 2 - female.\n";
            exit(0);
        }        
    }

    int numCompositeLoci = Loci.GetNumberOfCompositeLoci();

    vector<bool> true_vector(numCompositeLoci,true);
    _xi.assign(2,true_vector);
    sumxi.SetNumberOfElements(numCompositeLoci);

    vector<unsigned int> empty( 0, 0 );
    new_genotype.resize( numCompositeLoci, empty );
    PossibleHaplotypes = new Vector_i[numCompositeLoci];
    LocusAncestry.SetNumberOfElements( Loci.GetNumberOfChromosomes() );
    Matrix_d tm(1,1);
    ExpectedAncestry.resize( numCompositeLoci,tm );

    unsigned int lociI = 0;
    X_posn = 9999;
    string s1("\"X\"");
    for( unsigned int j = 0; j < numChromosomes; j++ ){
        if( chrm[j]->GetLabel(0) != s1 ){
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
        if( options->getPopulations() == 1 )
            LocusAncestry(j).SetElements(0);
        for( unsigned int jj = 0; jj < numCompLoci[j]; jj++ ){
            int compLocus = chrm[j]->GetLocus(jj);
            int numLoci = Loci(compLocus)->GetNumberOfLoci();
            vector<unsigned int> decodedGenotype(numLoci * 2,0);
            new_genotype[compLocus].resize( 2 * numLoci, 0 );
            for (int locus=0; locus<numLoci; locus++) {
                int allele0, allele1;

                if (options->IsPedFile() == 1) {
                    allele0 = StringConvertor::toInt(data[1 + options->genotypesSexColumn() + 2*lociI]);
                    allele1 = StringConvertor::toInt(data[2 + options->genotypesSexColumn() + 2*lociI]);
                } 
                else {
                    pair<int, int> a = StringConvertor::toIntPair(data[1 + options->genotypesSexColumn() + lociI]);
                    allele0 = a.first;
                    allele1 = a.second;
                }

                decodedGenotype[locus*2]   = allele0;
                decodedGenotype[locus*2+1] = allele1;
                new_genotype[compLocus][locus*2]      = allele0;
                new_genotype[compLocus][locus*2+1]    = allele1;
                lociI++;
            }
            _genotype.push_back(encodeGenotype(decodedGenotype));
        }
    }
    //may be possible to do this inside above loop but a loop over composite loci is neater
    for(int j=0;j<numCompositeLoci;++j)PossibleHaplotypes[j] = Loci(j)->SetPossibleHaplotypes(new_genotype[j]);
}

double
Individual::getLogLikelihoodXOnly( AdmixOptions* options, AlleleFreqs *A, Chromosome **chrm, Matrix_d ancestry, vector<double> rho )
{
   double LogLikelihood = 0.0;
   _rhoHat = rho;
   _ancestryHat = ancestry;
   Vector_d ff( A->GetNumberOfCompositeLoci() );
   vector< Vector_d > f(1,ff);
   for( unsigned int jj = 1; jj < numCompLoci[0]; jj++ ){
     f[0](jj) = exp( -A->getLoci()->GetDistance( jj ) * _rhoHat[0] );
   }
   chrm[0]->UpdateParametersHaploid( this, A,_ancestryHat, options, f, true );
   LogLikelihood += chrm[0]->getLogLikelihood();
   return LogLikelihood;
}
   
double
Individual::getLogLikelihood( AdmixOptions* options, AlleleFreqs *A, Chromosome **chrm, Matrix_d ancestry, vector<double> rho, Matrix_d ancestry_X, vector<double> rho_X )
{
   int locus = 0;
   _rhoHat = rho;
   _ancestryHat = ancestry;
   double LogLikelihood = 0.0;
   Vector_d ff( A->GetNumberOfCompositeLoci() );
   vector< Vector_d > f(2,ff);
   for( unsigned int j = 0; j < numChromosomes; j++ ){      
      locus++;
      if( j != X_posn ){
         for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
	   f[0](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat[0] );
            if( options->getModelIndicator() ){
	      f[1](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat[1] );
            }
            else
               f[1](locus) = f[0](locus);
            locus++;
         }
         chrm[j]->UpdateParameters( this,A, _ancestryHat, options, f, true);
         LogLikelihood += chrm[j]->getLogLikelihood();
      }
      else{
         _rhoHat_X = rho_X;
         _ancestryHat_X = ancestry_X;
         for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
	   f[0](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat_X[0] );
            if( sex == 2 ){
	      f[1](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rhoHat_X[1] );
            }
            locus++;
         }
         if( sex == 1 )
	   chrm[j]->UpdateParametersHaploid( this,A, _ancestryHat_X, options, f, true);
         else
	   chrm[j]->UpdateParameters( this, A,_ancestryHat_X, options, f, true);
         LogLikelihood += chrm[j]->getLogLikelihood();
      }
   }
   return LogLikelihood;
}

double
Individual::getLogLikelihoodOnePop(AlleleFreqs *A )
{
   double Likelihood = 0.0;
   Matrix_d Prob;
   for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
     if( _genotype[j][0] != 0 ){
       //CompositeLocus *locus = (CompositeLocus*)Loci(j);
       Prob = A->GetLikelihood( j, new_genotype[j], getPossibleHaplotypes(j), true, true );
        Likelihood += log( Prob(0,0) );
     }
   }
   return Likelihood;
}

double
Individual::getLogLikelihood( AdmixOptions* options, AlleleFreqs* A, Chromosome **chrm )
{
   int locus = 0;
   double LogLikelihood = 0.0;
   Vector_d ff( A->GetNumberOfCompositeLoci() );
   vector< Vector_d > f(2,ff);
   for( unsigned int j = 0; j < numChromosomes; j++ ){      
      locus++;
      if( j != X_posn ){
	for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
	  f[0](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho[0] );
            if( options->getModelIndicator() ){
	      f[1](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho[1] );
            }
            else
               f[1](locus) = f[0](locus);
            locus++;
         }
	chrm[j]->UpdateParameters( this, A,_ancestry, options, f, false);
         LogLikelihood += chrm[j]->getLogLikelihood();
      }
      else if( options->getXOnlyAnalysis() ){
         for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
	   f[0](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho[0] );
            locus++;
         }
         chrm[j]->UpdateParametersHaploid( this, A,_ancestry, options, f, false);
         LogLikelihood += chrm[j]->getLogLikelihood();
      }
      else{
         for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
	   f[0](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho_X[0] );
            if( sex == 2 ){
	      f[1](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho_X[1] );
            }
            locus++;
         }
         if( sex == 1 )
	   chrm[j]->UpdateParametersHaploid( this, A,_ancestry, options, f, false );
         else
	   chrm[j]->UpdateParameters( this, A,_ancestry, options, f, false );
         LogLikelihood += chrm[j]->getLogLikelihood();
      }
   }
   return LogLikelihood;
}

vector<unsigned int> Individual::encodeGenotype(vector<unsigned int>& decoded)
{
  int NumberOfLoci = decoded.size()/2;
  vector<unsigned int> encoded(decoded.size(),0);
  for(int i=0;i<NumberOfLoci;i++){
    encoded[0] *= 10;
    encoded[1] *= 10;
    encoded[0] += decoded[i*2];
    encoded[1] += decoded[i*2+1];
  }
  return encoded;
}

void
Individual::s2c(char *c, string s)
{
  int len = s.length();
  s.copy( c, len, 0 );
  c[len] = 0;
}


Individual::~Individual()
{
  delete []PossibleHaplotypes;

}

vector< unsigned int >& Individual::getGenotype(unsigned int locus)
{
  return new_genotype[locus];
}

Vector_i Individual::getPossibleHaplotypes(unsigned int locus){
  return PossibleHaplotypes[locus];
}

vector< vector< unsigned int > > & Individual::IsMissing()
{
  return _genotype;
}

vector< unsigned int >& Individual::IsMissing(unsigned int locus)
{
  return _genotype[locus];
}

void Individual::setGenotype(unsigned int locus,vector<unsigned int> genotype)
{//not used
  _genotype[locus] = genotype;
}

Matrix_d& Individual::getAncestry()
{
  //returns admixture proportions, should be renamed
  return _ancestry;
}

void Individual::setAncestry(Matrix_d ancestry)
{
  _ancestry = ancestry;
}

Matrix_d& Individual::getAncestryX()
{
  return _ancestryX;
}

void Individual::setAncestryX(Matrix_d ancestry)
{
  _ancestryX = ancestry;
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

Matrix_d Individual::getExpectedAncestry( int locus )
{
   return( ExpectedAncestry[ locus ] );
}

double Individual::getLogPosteriorProb()
{
   return LogPosterior;
}

// Samples individual admixture proportions conditional on sampled values of ancestry at loci where 
// jump indicator xi is 1, population admixture distribution parameters alpha, and likelihood from regression 
// model (if there is one)  
// should have an alternative function to sample population mixture component membership and individual admixture proportions
// conditional on genotypes, not sampled locus ancestry
Matrix_d
Individual::SampleParameters( int ind, AdmixOptions* options, AlleleFreqs *A, Chromosome **chrm, vector<Vector_d> alpha, bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta, int iteration, vector<double> sigma, Matrix_d &Theta_X )
{
   unsigned int locus = 0;
   sumxi.SetElements( 0 );
   //double q, Prob;
   Vector_d vectemp;
   vector< Vector_d > f( 2, A->getLociCorrSummary() );
   Matrix_i SumLocusAncestry( options->getPopulations(), 2 ), SumLocusAncestry_X;
   Matrix_d Theta;

   if( options->getModelIndicator() ){
      Theta.SetNumberOfElements( options->getPopulations(), 2 );
   }
   else{
      Theta.SetNumberOfElements( options->getPopulations(), 1 );
   }

   if( A->getLoci()->isX_data() ){
      if( sex == 1 ){
         Theta_X.SetNumberOfElements( options->getPopulations(), 1 );
         SumLocusAncestry_X.SetNumberOfElements( options->getPopulations(), 1 );
      }
      else{
         Theta_X.SetNumberOfElements( options->getPopulations(), 2 );
         SumLocusAncestry_X.SetNumberOfElements( options->getPopulations(), 2 );
      }
   }
  
  Sumrho0 = 0;
  // f0 and f1 are arrays of scalars of the form exp - rho*x, where x is distance between loci
  // required to calculate transition matrices 
  if( options->getRhoIndicator() ){
     for( unsigned int j = 0; j < numChromosomes; j++ ){
        locus++;
        for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
	  f[0](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho[0] );
           if( options->getModelIndicator() ){
	     f[1](locus) = exp( -A->getLoci()->GetDistance( locus ) * _rho[1] );
           }
           else
              f[1](locus) = f[0](locus);
           locus++;
        }
     }
  }

  SampleLocusAncestry(chrm, options, A, f, &SumLocusAncestry, &SumLocusAncestry_X);
  
  double L = A->getLoci()->GetLengthOfGenome(), L_X=0.0;
  if( A->getLoci()->isX_data() ) L_X = A->getLoci()->GetLengthOfXchrm();  
  vector< unsigned int > SumN(2,0);
  vector< unsigned int > SumN_X(2,0);
  if( options->getRhoIndicator() ){
 
    SampleNumberOfArrivals(A, options, &SumN, &SumN_X);

    SampleRho( options->getXOnlyAnalysis(), options->getModelIndicator(), A->getLoci()->isX_data(), rhoalpha, rhobeta, L, L_X, 
	       SumN, SumN_X);
  }
  
  if( options->getXOnlyAnalysis() ){
    vectemp = gendirichlet( alpha[0] + SumLocusAncestry_X.GetColumn(0) );
    Theta.SetColumn( 0, vectemp );
  }
  else if( options->getModelIndicator() ){
    for( unsigned int g = 0; g < 2; g++ ){
      if( options->getAnalysisTypeIndicator() > -1 )
	vectemp = gendirichlet( alpha[0] + SumLocusAncestry.GetColumn(g) );
      else
	vectemp = gendirichlet( alpha[g] + SumLocusAncestry.GetColumn(g) );
      Theta.SetColumn( g, vectemp );
    }
    if( A->getLoci()->isX_data() ){
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	vectemp = gendirichlet( Theta.GetColumn(g)*sigma[g]
				+ SumLocusAncestry_X.GetColumn(g) );
	Theta_X.SetColumn( g, vectemp );
      }
    }
  }
  else{
     vectemp = gendirichlet( alpha[0] + SumLocusAncestry.RowSum() );
     Theta.SetColumn( 0, vectemp );
  }
  if( options->getMLIndicator() && ind == 0 && iteration > options->getBurnIn() )
    CalculateLogPosterior(options,A->getLoci()->isX_data(), alpha, _symmetric,
				    _admixed,rhoalpha, rhobeta, L, L_X, 
				    SumN, SumN_X, SumLocusAncestry, SumLocusAncestry_X);

  return Theta;
}

void Individual::SampleLocusAncestry(Chromosome **chrm, AdmixOptions *options, AlleleFreqs *A, vector< Vector_d > f, 
				     Matrix_i *SumLocusAncestry, Matrix_i *SumLocusAncestry_X){
  // Loops over loci to sample locus ancestry and jump indicators xi 
  // computes conditional distribution of locus ancestry if required for score tests
  int locus=0;
  double q, Prob;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    locus = chrm[j]->GetLocus(0);

    if( j != X_posn ){
      chrm[j]->UpdateParameters( this, A,_ancestry, options, f, false );
      LocusAncestry(j) = chrm[j]->SampleForLocusAncestry( this, A );
    }
    else if( options->getXOnlyAnalysis() ){
      chrm[j]->UpdateParametersHaploid( this, A,_ancestry, options, f, false );
      LocusAncestry(j).SetRow( 0, chrm[j]->SampleForHaploidLocusAncestry( this, A ) );
    }
    else if( sex == 1 ){
      chrm[j]->UpdateParametersHaploid( this, A,_ancestryX, options, f, false );
      LocusAncestry(j).SetRow( 0, chrm[j]->SampleForHaploidLocusAncestry( this , A) );
    }
    else{
      chrm[j]->UpdateParameters( this, A,_ancestryX, options, f, false );
      LocusAncestry(j) = chrm[j]->SampleForLocusAncestry( this , A);
    }
    if( options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry())
      ExpectedAncestry[locus] = chrm[j]->getExpectedAncestry( 0 );
    locus++;
    for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
      if( options->getTestForAffectedsOnly()|| options->getTestForLinkageWithAncestry() )
	ExpectedAncestry[locus] = chrm[j]->getExpectedAncestry( jj );
      for( unsigned int g = 0; g < gametes[j]; g++ ){
	if( LocusAncestry(j)(g,jj-1) == LocusAncestry(j)(g,jj) ){
	  if( options->getModelIndicator() && g == 1 ){
	    q = _ancestry( LocusAncestry(j)(g,jj), g );
	  } else {
	    q = _ancestry( LocusAncestry(j)(g,jj), 0 );
	  }
	  Prob = (1 - f[g](locus)) / ( 1 - f[g](locus) + f[g](locus) / q );
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
      locus++;
    }
    locus = chrm[j]->GetLocus(0);
    // sum ancestry states over loci where jump indicator is 1
    for( unsigned int jj = 0; jj < numCompLoci[j]; jj++ ){
      for( unsigned int g = 0; g < gametes[j]; g++ ){
	if( _xi[g][locus] ){
	  if( j != X_posn )
	    (*SumLocusAncestry)( LocusAncestry(j)( g, jj ), g )++;
	  else
	    (*SumLocusAncestry_X)( LocusAncestry(j)( g, jj ), g )++;
	}
      }
      locus++;
    }
  }
}

void Individual::SampleNumberOfArrivals(AlleleFreqs *A, AdmixOptions *options, vector< unsigned int > *SumN, 
					vector< unsigned int > *SumN_X){
  // samples number SumN of arrivals between each pair of adjacent loci, 
  // conditional on jump indicators xi and sum of intensities rho
  // total number SumN is used for conjugate update of sum of intensities 
  double q, rho = 0.0;
  int locus = 0;
  int ran = 0;
  if( myrand() < 0.5 ) ran = 1;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    locus++;
    for( unsigned int jj = 1; jj < numCompLoci[j]; jj++ ){
      double delta = A->getLoci()->GetDistance(locus);
      for( unsigned int g = 0; g < gametes[j]; g++ ){
	if( _xi[g][locus] ){
	  if( options->getXOnlyAnalysis() ){
	    rho = _rho[0];
	    q = _ancestry( LocusAncestry(j)(0,jj-ran), 0 );
	  }
	  else if( options->getModelIndicator() && g == 1 ){
	    q = _ancestry( LocusAncestry(j)(g,jj-ran), g );
                    if( j != X_posn )
		      rho = _rho[g];
                    else
		      rho = _rho_X[g];
	  } else {
	    q = _ancestry( LocusAncestry(j)(g,jj-ran), 0 );
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
	    (*SumN)[g] += sample + 1;
	  else
	    (*SumN_X)[g] += sample + 1;
	}
      }
      locus++;
    }
  }
}

void Individual::SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta, double L, double L_X, 
			   vector< unsigned int > SumN, vector< unsigned int > SumN_X){
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

void Individual::CalculateLogPosterior(AdmixOptions *options, bool isX_data, vector<Vector_d> alpha, 
						 bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta, double L, 
						 double L_X, vector< unsigned int > SumN, vector< unsigned int > SumN_X, 
						 Matrix_i &SumLocusAncestry, Matrix_i &SumLocusAncestry_X){

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
      LogPosterior += getDirichletLogDensity(alphaparams0,_ancestryHat.GetColumn(0));
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
	LogPosterior += getDirichletLogDensity(alphaparams0,_ancestryHat.GetColumn(g));
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
	LogPosterior += getDirichletLogDensity(alphaparams0,_ancestryHat_X.GetColumn(g));
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
      x[0] += getDirichletLogDensity(alphaparams0,_ancestryHat.GetColumn(0));
      x[1] += getDirichletLogDensity(alphaparams0,_ancestryHat.GetColumn(1));
      alphaparams1 = alpha[1] + SumLocusAncestry.GetColumn(1);
      x[0] += getDirichletLogDensity(alphaparams1,_ancestryHat.GetColumn(1));
      x[1] += getDirichletLogDensity(alphaparams1,_ancestryHat.GetColumn(0));
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
	PR(_ancestryHat.GetColumn(0));
           PR(_ancestryHat.GetColumn(1));
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
           LogPosterior+=getDirichletLogDensity(alphaparams0,_ancestryHat.GetColumn(0));
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
	LogPosterior+=getDirichletLogDensity(alphaparams1,_ancestryHat.GetColumn(1));
      }
    }
  }
}

double Individual::IntegratingConst( double alpha, double beta, double a, double b )
{
   double I = gsl_cdf_gamma_P( b*beta, alpha, 1 ) - gsl_cdf_gamma_P( a*beta, alpha, 1);
   return I;
}

void Individual::UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates, Vector_d &poptheta, bool ModelIndicator,
Matrix_d *Covariates0)
{
  Vector_d avgtheta;
  if(ModelIndicator )
    avgtheta = _ancestry.RowMean();
  else
    avgtheta = _ancestry.GetColumn(0);
  for( int k = 0; k < Populations - 1; k++ )
    (*Covariates0)( i, NoCovariates - Populations + k + 1 )
      = avgtheta( k + 1 ) - poptheta( k + 1 );
}


void Individual::Accept_Reject_Theta( double p, Matrix_d &theta, Matrix_d &thetaX, bool xdata, int Populations, bool ModelIndicator )
{
  bool test = true;
  // loop over populations: if element of Dirichlet parameter vector is 0, do not update corresponding element of 
  // admixture proportion vector
  for( int k = 0; k < Populations; k++ ){
    if( (theta)( k, 0 ) == 0.0 )
      test = false;
    else if( ModelIndicator && (theta)( k, 1 ) == 0.0 )
      test = false;
  }

  // generic Metropolis rejection step
  if( p < 0 ){
     if( log(myrand()) < p && test ){
        setAncestry(theta);
        if( xdata )
           setAncestryX(thetaX);
     }
  }
  else{
     setAncestry(theta);
     if( xdata )
        setAncestryX(thetaX);
  }
}

// these two functions return ratio of likelihoods of new and old values of population admixture
// in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in. 
 
// should have just one function to get the likelihood in the regression model, given a value of population admixture
// should be generic method for GLM, given Xbeta, Y and probability distribution 
// then should calculate ratio in Metropolis step 
//  
double Individual::AcceptanceProbForTheta_LogReg( int i, int TI, Matrix_d &theta ,bool ModelIndicator,int Populations,
					      int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, MatrixArray_d &ExpectedY, 
						  MatrixArray_d &Target, Vector_d &poptheta) 
{
  double prob, Xbeta = 0;
  // TI appears to define which outcome var is used
  Vector_d avgtheta;
  // calculate mean of parental admixture proportions
  if( ModelIndicator )
    avgtheta = theta.RowMean() - poptheta;
  else
    avgtheta = theta.GetColumn(0) - poptheta;

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

double Individual::AcceptanceProbForTheta_LinearReg( int i, int TI,  Matrix_d &theta ,bool ModelIndicator,int Populations,
						 int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, MatrixArray_d &ExpectedY,
						 MatrixArray_d &Target, Vector_d &poptheta, Vector_d &lambda)
{
  double prob, Xbeta = 0;
  Vector_d avgtheta;
  if( ModelIndicator )
    avgtheta = theta.RowMean() - poptheta;
  else
    avgtheta = theta.GetColumn(0) - poptheta;

  for( int jj = 0; jj < NoCovariates - Populations + 1; jj++ )
    Xbeta += Covariates0( i, jj ) * beta( TI )( jj, 0 );
  for( int k = 0; k < Populations - 1; k++ ){
    Xbeta += avgtheta( k ) * beta( TI )( NoCovariates - Populations + k + 1, 0 );
  }

  prob = 0.5 * lambda( TI ) * (( ExpectedY( TI )( i, 0 ) - Target( TI )( i, 0 ) ) * ( ExpectedY( TI )( i, 0 ) - Target( TI )( i, 0 ) )
			       - ( Xbeta - Target( TI )( i, 0 ) ) * ( Xbeta - Target( TI )( i, 0 ) ) );

  return( prob );
}

double Individual::AcceptanceProbForTheta_XChrm(Matrix_d &Theta, Matrix_d &ThetaX,std::vector<double> &sigma, int Populations )
{
   int gametes = 1;
   if( sex == 2 )
      gametes = 2;
   double p = 0;
   Matrix_d ThetaOld = _ancestry, ThetaXOld = _ancestryX;
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
					     double rhobeta,vector<double> sigma)
{
  double u;
  Matrix_d Theta, ThetaX;

  if( options->getAnalysisTypeIndicator() > 1 ){
    for( int k = 0; k < Target->GetNumberOfElements(); k++ ){
      if( (*Target)(k).IsMissingValue( i, 0 ) ){
	if( !OutcomeType(k) )
	  (*Target)(k)( i, 0 ) = gennor( ExpectedY(k)( i, 0 ), 1 / sqrt( lambda(k) ) );
	else{
	  u = myrand();
	  if( u * ExpectedY(k)( i, 0 ) < 1 )
	    (*Target)(k)( i, 0 ) = 1;
	  else
	    (*Target)(k)( i, 0 ) = 0;
	}
      }
    }
  }
  // should be modified to allow a population mixture component model   
  Theta = SampleParameters(i, options, A, chrm, alpha, _symmetric, _admixed,rhoalpha, rhobeta, iteration, sigma, ThetaX );

   double p = 0;
   if( options->getAnalysisTypeIndicator() == 2 && !options->getScoreTestIndicator() ){
     p = AcceptanceProbForTheta_LinearReg( i, 0, Theta ,options->getModelIndicator(),options->getPopulations(),
					      NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta,lambda); 
   }
   else if( (options->getAnalysisTypeIndicator() == 3 || options->getAnalysisTypeIndicator() == 4) && !options->getScoreTestIndicator() ){
     p = AcceptanceProbForTheta_LogReg( i, 0, Theta ,options->getModelIndicator(),options->getPopulations(),
						 NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta); 
   }
   else if( options->getAnalysisTypeIndicator() == 5 ){
      for( int k = 0; k < Target->GetNumberOfElements(); k++ ){
         if( OutcomeType( k ) )
	   p += AcceptanceProbForTheta_LogReg( i, k, Theta, options->getModelIndicator(), options->getPopulations(),
						 NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta); 
         else
	   p += AcceptanceProbForTheta_LinearReg( i, k, Theta, options->getModelIndicator(), options->getPopulations(),
					      NoCovariates, Covariates0, beta, ExpectedY, *Target, poptheta,lambda);
      }
   }
   if( A->getLoci()->isX_data() && !options->getXOnlyAnalysis() )
     p += AcceptanceProbForTheta_XChrm( Theta, ThetaX , sigma, options->getPopulations());
   Accept_Reject_Theta(p, Theta, ThetaX, A->getLoci()->isX_data(),options->getPopulations(), options->getModelIndicator() );
   if( options->getAnalysisTypeIndicator() > 1 )
     UpdateAdmixtureForRegression(i,options->getPopulations(), NoCovariates, poptheta, options->getModelIndicator(),&(Covariates0));
   for( int k = 0; k < options->getPopulations(); k++ ){
     (*SumLogTheta)( k ) += log( _ancestry( k, 0 ) );
     if(options->getModelIndicator() && !options->getXOnlyAnalysis() )
       (*SumLogTheta)( k ) += log( _ancestry( k, 1 ) );
   }
}

void Individual::OnePopulationUpdate( int i, MatrixArray_d *Target, Vector_i &OutcomeType, MatrixArray_d &ExpectedY, Vector_d &lambda, 
				     AlleleFreqs *A, int AnalysisTypeIndicator )
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
      
  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
    A->UpdateAlleleCounts(j,getGenotype(j),getPossibleHaplotypes(j), ancestry );
  }
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
			       getAncestry().GetColumn(0))
	+getDirichletLogDensity(alpha[1],
				getAncestry().GetColumn(1));
      if( *LogLikelihood > MaxLogLikelihood[i] ){
	*LogFileStreamPtr << getAncestry()
			 << _rho[0] << " " << _rho[1]
			 << endl << *LogLikelihood << endl
			 << iteration << endl;
	MaxLogLikelihood[i] = *LogLikelihood;
	if( iteration <= options->getBurnIn() ){
	  thetahat(i) = getAncestry();
	  rhohat[i] = _rho;
	  A->setAlleleFreqsMAP();
	  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
	    CompositeLocus *locus = (CompositeLocus*)(A->getLocus(j));
	    //locus->setAlleleFreqsMAP();
	    if( locus->GetNumberOfLoci() > 2 )
	      locus->setHaplotypeProbsMAP();
	  }
	  if( A->getLoci()->isX_data() ){
	    thetahatX(i) = getAncestryX();
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
