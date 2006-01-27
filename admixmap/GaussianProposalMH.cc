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
  newnum = *x;
  NewtonRaphson(args);//sets ddf, 2nd derivative at mode, used to set precision of proposal distribution

  //generate proposal from Normal distribution centred on mode of target distribution with precision 
  //given as sqrt of 2nd derivative at the mode
  double xnew = gennor( newnum, 1 / sqrt( -ddf ) );

  double LogTarget = (*function)( *x, args );//log posterior at current value
  double NewLogTarget = (*function)( xnew, args );//log posterior at proposal

  double ProposalRatio = LogNormalDensity( xnew, newnum, -ddf ) - LogNormalDensity( *x, newnum, -ddf );

  double LogAcceptanceProb = NewLogTarget - LogTarget + ProposalRatio;
  int flag = 0;
  if( log( myrand() ) < LogAcceptanceProb ){
    *x = xnew;
    flag = 1;
  }
  return flag;//returns indicator variable for acceptance
}

void GaussianProposalMH::NewtonRaphson(const void* const args)
{
  double step, df;
  do{
    ddf = (*ddfunction)( newnum, args );//2nd derivative
    df = (*dfunction)( newnum, args );//1st derivative
    step = -df / ddf;
    newnum += step;
  }while( fabs(df) > 0.001 );
}

double GaussianProposalMH::LogNormalDensity(double x, double mu, double lambda)
//returns log density of Normal distribution with mean mu and precision lambda at x
{
  return( -0.5 * lambda * ( x - mu ) * ( x - mu ) );
}
