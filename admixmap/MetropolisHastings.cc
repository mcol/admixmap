#include "MetropolisHastings.h"

MetropolisHastings::MetropolisHastings
(const Vector_d &inparameters,
 double (*funct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
 double (*dfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
 double (*ddfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
 const MatrixArray_i &integer_data, const MatrixArray_d &double_data )
{
   parameters = inparameters;
   data_i = integer_data;
   data_d = double_data;
   function = funct;
   dfunction = dfunct;
   ddfunction = ddfunct;
}

MetropolisHastings::~MetropolisHastings()
{
}

void MetropolisHastings::UpdateParameters( const Vector_d &inparameters )
{
  parameters = inparameters;
}

void MetropolisHastings::UpdateIntegerData( const MatrixArray_i &indata )
{
  data_i = indata;
}

void MetropolisHastings::UpdateDoubleData( const MatrixArray_d &indata )
{
  data_d = indata;
}

int MetropolisHastings::Sample( double *x )
{
  int flag = 0;
  double xnew, LogPost, NewLogPost ,ProposalRatio, LogAcceptanceProb;
  newnum = *x;
  NewtonRaphson();

  xnew = gennor( newnum, 1 / sqrt( -ddf ) );

  LogPost = (*function)( parameters, data_i, data_d, *x );
  NewLogPost = (*function)( parameters, data_i, data_d, xnew );

  ProposalRatio = LogNormalDensity( xnew, newnum, -ddf ) - LogNormalDensity( *x, newnum, -ddf );

  LogAcceptanceProb = NewLogPost - LogPost + ProposalRatio;
  if( log( myrand() ) < LogAcceptanceProb ){
    *x = xnew;
    flag = 1;
  }
  return flag;
}

void MetropolisHastings::NewtonRaphson()
{
  double step, df;
  do{
    ddf = (*ddfunction)( parameters, data_i, data_d, newnum );
    df = (*dfunction)( parameters, data_i, data_d, newnum );
    step = -df / ddf;
    newnum += step;
  }while( fabs(df) > 0.001 );
}

double MetropolisHastings::LogNormalDensity
(double x, double mu, double lambda)
{
  return( -0.5 * lambda * ( x - mu ) * ( x - mu ) );
}
