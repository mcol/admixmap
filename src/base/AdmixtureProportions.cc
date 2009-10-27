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

#include "AdmixtureProportions.h"
#include "PedBase.h"

using bclib::pvector;


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
