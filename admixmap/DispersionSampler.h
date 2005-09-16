// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   DispersionSampler.h 
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

class DispersionSampler{

public:
  DispersionSampler();
  ~DispersionSampler();

  void setDimensions(unsigned inL, unsigned inK, int* const inH, double step0, double min, double max, double target);
  void setEtaPrior(double, double);
  void addAlphas(unsigned, const double* const);
  void addCounts(unsigned, const int* const);
  double Sample();
  double getEnergy(double);
  double getGradient(double);

  float getAcceptanceRate();
  float getStepsize();

private:
  //dimensions
  unsigned L;
  unsigned K;//number of counts

  double logeta[1];
  HamiltonianMonteCarlo Sampler;
  double **Args;

  static double etaEnergyFunction(unsigned , const double * const logitmu, const double* const *args);
  static void etaGradient(unsigned , const double * const logitmu, const double* const *args, double *g);
};
