/** 
 *   ADMIXMAP
 *   GaussianProposalMH.cc (formerly MetropolisHastings.cc)
 *   This class is used to implement a Metropolis Hastings update with Gaussian proposal distribution
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
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
  //determine mode of target distribution
  //and set ddf to 2nd derivative at mode
  double newnum = FindModeByNewtonRaphson(*x, args);

  //generate proposal from Normal distribution centred on mode of target distribution with precision 
  //given as sqrt of -2nd derivative at the mode
  double xnew = Rand::gennor( newnum, 1 / sqrt( -ddf ) );

  double LogTarget = (*function)( *x, args );//log of target dist at current value
  double NewLogTarget = (*function)( xnew, args );//log of target dist at proposal

  double ProposalRatio = LogNormalDensity( xnew, newnum, -ddf ) - LogNormalDensity( *x, newnum, -ddf );

  double LogAcceptanceProb = NewLogTarget - LogTarget + ProposalRatio;
  int accept = 0;
  if( log( Rand::myrand() ) < LogAcceptanceProb ){
    *x = xnew;
    accept = 1;
  }

  return accept;//returns indicator variable for acceptance
}

double GaussianProposalMH::FindModeByNewtonRaphson(double init, const void* const args)
{
  //init is an initial guess, args are arguments to loglikelihood and derivative functions
  //returns mode as estimate of mode of f
  double step, f, flast;
  
  double df = (*dfunction)( init, args );//1st derivative at initial value
  ddf = (*ddfunction)( init, args );//2nd derivative " " 

  double newnum = init;
  while( fabs(df) > 0.001 ){
    flast = (*function)(newnum, args);//loglikelihood at current value
    step = -df / ddf;

    f = (*function)(step+newnum, args);//loglikelihood at new value
    if(f > flast)//if NR step increases loglikelihood
      newnum += step;//use it
    else{//use simple random walk to find a point nearer the mode
      //flast = f; 
      do{
	step = Rand::gennor(newnum, 0.1 / sqrt( -ddf ) );//could improve this, maybe with StepSizeTuner
	f = (*function)(step, args);
      }
      while(f < flast);//repeat until loglikelihood increases
      newnum = step;
    }

    df = (*dfunction)( newnum, args );//update 1st derivative
    ddf = (*ddfunction)( newnum, args );//update 2nd derivative
  }
  return newnum;
}

double GaussianProposalMH::LogNormalDensity(double x, double mu, double lambda)
//returns log density of Normal distribution with mean mu and precision lambda at x
{
  return( -0.5 * lambda * ( x - mu ) * ( x - mu ) );
}
