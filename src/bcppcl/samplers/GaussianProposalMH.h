// *-*-C++-*-*
/** 
 *   GaussianProposalMH.h
 *   This class is used to implement a Metropolis Hastings update with Gaussian proposal distribution
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef GAUSSIANPROPOSALMH_H
#define GAUSSIANPROPOSALMH_H 1

///Class to implement a Metropolis Hastings update with Gaussian proposal distribution
class GaussianProposalMH
{
public:
  GaussianProposalMH(double (*funct)( const double, const void* const ),
                     double (*dfunct)( const double, const void* const ),
                     double (*ddfunct)( const double, const void* const )
                     );
   
  ~GaussianProposalMH();
  
  int Sample(double*, const void* const);

private: // members
  unsigned dim;
  double ddf;

  double (*function)( const double, const void* const );
  double (*dfunction)( const double, const void* const );
  double (*ddfunction)( const double, const void* const );

  // HELPER FUNCTIONS
  double FindModeByNewtonRaphson(double init, const void* const args);
  static double LogNormalDensity
  (double x, double mu, double lambda);


  // NOT IMPLEMENTED!!
  // - to prevent use
  // default constructor
  GaussianProposalMH();
  // private copy constructor
  GaussianProposalMH(const GaussianProposalMH&);
  // private assignment operator
  GaussianProposalMH& operator=(const GaussianProposalMH&);

};

#endif /* !defined GAUSSIANPROPOSALMH_H */
