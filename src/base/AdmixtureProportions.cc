/*
 *  Copyright (c) 2009 Marco Colombo
 *
 *  This is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this software; see the file COPYING. If not, it can be
 *  found at http://www.gnu.org/copyleft/gpl.html or by writing to the
 *  Free Software Foundation, 51 Franklin Street, Fifth Floor, Boston,
 *  MA 02110-1301, USA.
 */

/**
 *  @file AdmixtureProportions.cc
 *  Implementation of the AdmixtureProportions class.
 */

#include <cmath>
#include "AdmixtureProportions.h"
#include "PedBase.h"

using genepi::cvector;
using bclib::pvector;

static void oddsratios2xKtocolratios(const cvector<double>& psi,
                                     const pvector<double>& P,
                                     double Q, cvector<double>& r);

/// Constructor
AdmixtureProportions::AdmixtureProportions(int _nGametes, int _nHiddenStates) :
  gametes(_nGametes), vflat(NULL), isOk(true) {

  gp_assert(_nGametes >= 1);
  gp_assert(_nHiddenStates >= 1);

  for (int i = 0; i < _nGametes; ++i)
    gametes[i].resize(_nHiddenStates);
}

/// Copy constructor
AdmixtureProportions::AdmixtureProportions(const AdmixtureProportions& rhs) :
  gametes(rhs.gametes.size()), vflat(NULL), isOk(true) {

  const size_t G = gametes.size();
  const size_t K = rhs.gametes[0].size();
  for (size_t i = 0; i < G; ++i)
    gametes[i].resize(K);

  for (size_t g = 0; g < G; ++g) {
    for (size_t k = 0; k < K; ++k)
      gametes[g][k] = rhs[g][k];
  }
}

/// Destructor
AdmixtureProportions::~AdmixtureProportions() {

  delete[] vflat;
}

/// Assignment operator
AdmixtureProportions& AdmixtureProportions::operator=(const AdmixtureProportions& rhs) {

  if (this != &rhs) {

    const size_t G = gametes.size();
    const size_t K = gametes[0].size();
    gp_assert(rhs.gametes.size() == G);
    gp_assert(rhs[0].size() == K);

    for (size_t g = 0; g < G; ++g) {
      for (size_t k = 0; k < K; ++k)
        gametes[g][k] = rhs[g][k];
    }
  }

  return *this;
}

/// Set the dimensions of an existing object
void AdmixtureProportions::setDimensions(int _nGametes, int _nHiddenStates) {

  gp_assert(_nGametes >= 1);
  gp_assert(_nHiddenStates >= 1);

  gametes.resize(_nGametes);
  for (int i = 0; i < _nGametes; ++i)
    gametes[i].resize(_nHiddenStates);

  if (vflat) {
    delete[] vflat;
    vflat = NULL;
  }

  isOk = true;
}

/// Get direct write access to the pvector associated to the requested gamete
pvector<double>& AdmixtureProportions::operator[](int idx) {

  gp_assert(isOk);
  return gametes[idx];
}

/// Get direct read access to the pvector associated to the requested gamete
const pvector<double>& AdmixtureProportions::operator[](int idx) const {

  gp_assert(isOk);
  return gametes[idx];
}

/// Set all elements to the specified value
void AdmixtureProportions::setTo(double value) {

  gp_assert(isOk);

  const size_t G = gametes.size();
  const size_t K = gametes[0].size();

  for (size_t g = 0 ; g < G; ++g)
    for (size_t k = 0 ; k < K ; ++k)
      gametes[g][k] = value;
}

/// Return a (read-only) flattened version of the internal vectors
const double* AdmixtureProportions::flat() const {

  gp_assert(isOk);

  const size_t G = gametes.size();
  const size_t K = gametes[0].size();

  if (!vflat)
    vflat = new double[G * K];

  for (size_t g = 0 ; g < G; ++g) {
    for (size_t k = 0 ; k < K ; ++k)
      vflat[g * K + k] = gametes[g][k];
  }

  return vflat;
}

/// Return a (read-only) flattened version of the internal vectors adjusted
/// to account for the admixture ratio of the X chromosome
const double* AdmixtureProportions::flatXChromosome(const cvector<double>& psi) const {

  gp_assert(isOk);
  gp_assert(psi.size() == gametes[0].size());

  const size_t G = gametes.size();
  const size_t K = gametes[0].size();
  const double Q = 0.5;
  cvector<double> r(K);
  double value;

  if (!vflat)
    vflat = new double[G * K];

  for (size_t g = 0 ; g < G; ++g) {

    oddsratios2xKtocolratios(psi, gametes[g], Q, r);

    for (size_t k = 0 ; k < K ; ++k) {
      value = 2.0 * gametes[g][k] * (1.0 + 2.0 * r[k]) / (3.0 * (1.0 + r[k]));
      vflat[g * K + k] = value;
    }
  }

  return vflat;
}

/// Print the content of the internal vectors
void AdmixtureProportions::print() const {

  gp_assert(isOk);

  const size_t G = gametes.size();
  const size_t K = gametes[0].size();

  for (size_t g = 0 ; g < G; ++g) {
    for (size_t k = 0 ; k < K ; ++k)
      printf("Theta[%zd][%zd] = %f\n", g, k, gametes[g][k]);
  }
}

/// Given the odds ratios psi and a vector of frequencies P, compute the
/// ratio of paternal founders to maternal founders for each population
void oddsratios2xKtocolratios(const cvector<double>& psi,
                              const pvector<double>& P,
                              double Q, cvector<double>& r) {

  // This implements a regula falsi bisection algorithm following
  //    http://en.wikipedia.org/wiki/False_position_method
  //
  // The equation to solve is
  //    f(x) = -Q + sum_i p_i / (1 + psi_i * exp(x)) = 0
  // which has the convenient property that
  //    lim x->+inf = -Q       < 0
  //    lim x->-inf = sum_i pi > 0
  // so that by choosing a sufficently large value of x, we can produce
  // an initial bracketing of the root to initialise the bisection.
  // The value must be large enough to guarantee that each endpoint of
  // the bracketing interval has opposite sign; however, too large values
  // may slow down the convergence to the root.

#ifdef DEBUG_BISECTION
  int iter = 0;
  for (size_t i = 0; i < P.size(); ++i)
    gp_assert(psi[i] > 0.0);
#endif

  gp_assert(Q > 0.0);
  gp_assert(Q < 1.0);

  const int K = (int) P.size();

  // parameters for the bisection
  const double epsilon = 1.0e-9;
  const double limit = 15.0;
  double left = -limit, right = limit, diff = right - left;

  // store which side of the interval has been updated
  int side = 0;

  // f(x) = sum_i ( p_i / (1 + psi_i * exp(x)) ) - Q
  double fleft  = -Q, expleft  = exp(left);
  double fright = -Q, expright = exp(right);
  for (int i = 0; i < K; ++i) {
    fleft  += P[i] / (1.0 + psi[i] * expleft);
    fright += P[i] / (1.0 + psi[i] * expright);
  }

  double mid, fmid, expmid;

  while (fabs(diff) > epsilon) {

    // evaluate the function at the mid point found through the regula falsi
    mid = (left * fright - right * fleft) / (fright - fleft);
    fmid = -Q, expmid = exp(mid);
    for (int i = 0; i < K; ++i)
      fmid += P[i] / (1.0 + psi[i] * expmid);

#ifdef DEBUG_BISECTION
    printf("it %d: bisect on [%g, %g] -> %g (diff: %.2e)\n",
           iter, left, right, mid, diff);
#endif

    // update the endpoints of the bisection interval
    if ((fleft > 0 && fmid > 0) || (fleft < 0 && fmid < 0)) {
      left = mid;
      fleft = fmid;
      if (side == 1)
        fright = fright * 0.5;
      side = 1;
    }
    else if ((fright > 0 && fmid > 0) || (fright < 0 && fmid < 0)) {
      right = mid;
      fright = fmid;
      if (side == -1)
        fleft = fleft * 0.5;
      side = -1;
    }
    else {

#ifdef DEBUG_BISECTION
      if (fabs(fmid) > epsilon) {
        // this happens for values of psi that tend to infinity or to zero
        // the workaround is to enlarge the initial [left, right] interval
        printf("\n*** Bisection failed! fmid=%g\n", fmid);
        printf("Psi: ");
        for (int i = 0; i < K; ++i)
          printf("%f ", psi[i]);
        printf("\n");
        printf("P: ");
        for (int i = 0; i < K; ++i)
          printf("%f ", P[i]);
        printf("\n");
      }
#endif
      gp_assert(fabs(fmid) < 2 * epsilon);
      break;
    }

    diff = right - left;

#ifdef DEBUG_BISECTION
    printf("   it %d: y(%f)=%e diff=%e\n", iter++, mid, fmid, diff);
#endif
  }

  expmid = exp(mid);
  for (int i = 0; i < K; ++i)
    r[i] = expmid * psi[i];
}
