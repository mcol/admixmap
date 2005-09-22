// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   MuSampler.h 
 *   Class to sample the proportion parameters of a multinomial-Dirichlet distribution
 *   parameterised as \mu_1, ..., \mu_H, \eta ie a vector of proportions and a dispersion parameter
 *   Copyright (c) 2005 LSHTM
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

#include "HamiltonianMonteCarlo.h"
#include "StepSizeTuner.h"
#include "DARS.h"

class MuSampler{

public:
  MuSampler();
  ~MuSampler();

  void setDimensions(unsigned inK, unsigned inH, double mustep0, double mumin, double mumax, double mutarget);
  void Sample(double* const alpha, double eta, const int* const Counts);
  void Sample1D(double* alpha, double eta, const int* const Counts);

  float getAcceptanceRate();
  float getStepsize();

private:
  //dimensions
  unsigned K;//number of counts
  unsigned H;//dimension of alpha

  //proportions are transformed using softmax transformation and sampled with a Hamiltonian Monte Carlo sampelr
  double *params;
  HamiltonianMonteCarlo muSampler;
  double **muArgs;

  static double muEnergyFunction(unsigned , const double * const params, const double* const *args);
  static void muGradient(unsigned , const double * const params, const double* const *args, double *g);
  static double fMu( const double* parameters, const int *counts,  const double *, double alpha );
  static double dfMu( const double* parameters, const int *counts, const double *, double alpha );
  static double ddfMu( const double* parameters, const int *counts, const double *, double alpha );
  static double logJacobian(const double* a, const double z, unsigned H);
  static double DlogJacobian(const double* const a, const double z, unsigned H, unsigned h, double delta);
};
