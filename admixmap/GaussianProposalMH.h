// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   GaussianProposalMH.h (formerly MetropolisHastings.h)
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
#ifndef GAUSSIANPROPOSALMH_H
#define GAUSSIANPROPOSALMH_H 1

#include "rand.h"
#include "matrix_i.h"
#include "matrix_d.h"

class GaussianProposalMH
{
public:
  GaussianProposalMH(const double* inparameters,
                     double (*funct)(const double*, Matrix_i&, Matrix_d&, double),
                     double (*dfunct)(const double*, Matrix_i&, Matrix_d&, double),
                     double (*ddfunct)(const double*, Matrix_i&, Matrix_d&, double),
                     const Matrix_i&, const Matrix_d&);
   
  ~GaussianProposalMH();
  
  int Sample(double*);
  void UpdateParameters(const double* inparameters);
  void UpdateIntegerData(const Matrix_i&);
  void UpdateDoubleData(const Matrix_d&);

private: // members
  unsigned dim;
  const double *parameters;
  Matrix_i data_i;
  Matrix_d data_d;
  double newnum;
  double ddf;
  double (*function)(const double*, Matrix_i &, Matrix_d &, double xin);
  double (*dfunction)(const double*, Matrix_i &, Matrix_d &, double xin);
  double (*ddfunction)(const double*, Matrix_i &, Matrix_d &, double xin);

  // HELPER FUNCTIONS
  void NewtonRaphson();
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
