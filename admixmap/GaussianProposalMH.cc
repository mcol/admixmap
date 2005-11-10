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
(double (*funct)( const double, const void* const ),
 double (*dfunct)( const double, const void* const ),
 double (*ddfunct)( const double, const void* const )
 ) 
{
   function = funct;
   dfunction = dfunct;
   ddfunction = ddfunct;
}

GaussianProposalMH::~GaussianProposalMH()
{
}

int GaussianProposalMH::Sample( double *x, const void* const args )
{
  int flag = 0;
  double xnew, LogPost, NewLogPost ,ProposalRatio, LogAcceptanceProb;
  newnum = *x;
  NewtonRaphson(args);//sets ddf

  xnew = gennor( newnum, 1 / sqrt( -ddf ) );

  LogPost = (*function)( *x, args );
  NewLogPost = (*function)( xnew, args );

  ProposalRatio = LogNormalDensity( xnew, newnum, -ddf ) - LogNormalDensity( *x, newnum, -ddf );

  LogAcceptanceProb = NewLogPost - LogPost + ProposalRatio;
  if( log( myrand() ) < LogAcceptanceProb ){
    *x = xnew;
    flag = 1;
  }
  return flag;
}

void GaussianProposalMH::NewtonRaphson(const void* const args)
{
  double step, df;
  do{
    ddf = (*ddfunction)( newnum, args );
    df = (*dfunction)( newnum, args );
    step = -df / ddf;
    newnum += step;
  }while( fabs(df) > 0.001 );
}

double GaussianProposalMH::LogNormalDensity(double x, double mu, double lambda)
{
  return( -0.5 * lambda * ( x - mu ) * ( x - mu ) );
}
