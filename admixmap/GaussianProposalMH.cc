/** 
 *   ADMIXMAP
 *   GaussianProposalMH.cc (formerly MetropolisHastings.cc)
 *   This class is used to implement a Metropolis Hastings update with Gaussian proposal distribution
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include "GaussianProposalMH.h"

GaussianProposalMH::GaussianProposalMH
(const double *inparameters,
 double (*funct)(const double*, Matrix_i&, Matrix_d&, double),
 double (*dfunct)(const double*, Matrix_i&, Matrix_d&, double),
 double (*ddfunct)(const double*, Matrix_i&, Matrix_d&, double),
 const Matrix_i &integer_data, const Matrix_d &double_data )
{
   parameters = inparameters;
   data_i = integer_data;
   data_d = double_data;
   function = funct;
   dfunction = dfunct;
   ddfunction = ddfunct;
}

GaussianProposalMH::~GaussianProposalMH()
{
}

void GaussianProposalMH::UpdateParameters( const double *inparameters )
{//may be unnecessary
  parameters = inparameters;
}

void GaussianProposalMH::UpdateIntegerData( const Matrix_i &indata )
{
  data_i = indata;
}

void GaussianProposalMH::UpdateDoubleData( const Matrix_d &indata )
{
  data_d = indata;
}

int GaussianProposalMH::Sample( double *x )
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

void GaussianProposalMH::NewtonRaphson()
{
  double step, df;
  do{
    ddf = (*ddfunction)( parameters, data_i, data_d, newnum );
    df = (*dfunction)( parameters, data_i, data_d, newnum );
    step = -df / ddf;
    newnum += step;
  }while( fabs(df) > 0.001 );
}

double GaussianProposalMH::LogNormalDensity
(double x, double mu, double lambda)
{
  return( -0.5 * lambda * ( x - mu ) * ( x - mu ) );
}
