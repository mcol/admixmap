/** 
 *   GaussianProposalMH.cc
 *   This class is used to implement a Metropolis Hastings update with Gaussian proposal distribution
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "bclib/GaussianProposalMH.h"
#include "bclib/ModeFinder.h"
#include <math.h>
#include "bclib/rand.h"

BEGIN_BCLIB_NAMESPACE

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
  //determine mode of target distribution
  //and set ddf to 2nd derivative at mode
  //double mode = FindModeByNewtonRaphson(*x, args);
  double mode = bclib::ModeFinder::FindModeByNewtonRaphson(*x, &ddf, args, function, dfunction, ddfunction);

  //generate proposal from Normal distribution centred on mode of target distribution with precision 
  //given as sqrt of -2nd derivative at the mode
  double xnew = bclib::Rand::gennor( mode, 1 / sqrt( -ddf ) );

  double LogTarget = (*function)( *x, args );//log of target dist at current value
  double NewLogTarget = (*function)( xnew, args );//log of target dist at proposal

  double ProposalRatio = LogNormalDensity( xnew, mode, -ddf ) - LogNormalDensity( *x, mode, -ddf );

  double LogAcceptanceProb = NewLogTarget - LogTarget + ProposalRatio;
  int accept = 0;
  if( log( bclib::Rand::myrand() ) < LogAcceptanceProb ){
    *x = xnew;
    accept = 1;
  }

  return accept;//returns indicator variable for acceptance
}

//TODO: move out
double GaussianProposalMH::LogNormalDensity(double x, double mu, double lambda)
//returns log density of Normal distribution with mean mu and precision lambda at x
{
  return( -0.5 * lambda * ( x - mu ) * ( x - mu ) );
}

END_BCLIB_NAMESPACE
