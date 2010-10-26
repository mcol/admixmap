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

//=============================================================================
/// \file AdmixtureProportions.h
/// Definition of the AdmixtureProportions class.
//=============================================================================

#ifndef ADMIXTUREPROPORTIONS_H
#define ADMIXTUREPROPORTIONS_H 1

#include <bclib/cvector.h>
#include <bclib/pvector.h>


/** \addtogroup base
 * @{ */


typedef genepi::cvector< bclib::pvector<double> > ThetaType;
typedef genepi::cvector< double >                 PsiType;


/// Admixture proportions
class AdmixtureProportions {

 public:

  /// Default constructor (empty)
  AdmixtureProportions() : vflat(NULL), isOk(false) { };

  /// Constructor
  AdmixtureProportions(int _nGametes, int _nPopulations);

  /// Copy constructor
  AdmixtureProportions(const AdmixtureProportions& rhs);

  /// Destructor
  ~AdmixtureProportions();

  /// Assignment operator
  AdmixtureProportions& operator=(const AdmixtureProportions& rhs);

  /// Set the dimensions of an existing object
  void setDimensions(int _nGametes, int _nPopulations);

  /// Return the number of gametes for which we have allocated space
  size_t getNGametes() const;

  /// Get direct write access to the pvector associated to the requested gamete
  bclib::pvector<double>& operator[](int idx);

  /// Get direct read access to the pvector associated to the requested gamete
  const bclib::pvector<double>& operator[](int idx) const;

  /// Set all elements to the specified value
  void setTo(double val);

  /// Scale all elements by the specified value
  void scaleBy(double val);

  /// Return a copy of the internal vectors
  ThetaType getTheta() const;

  /// Return a copy of the internal vectors, adjusted to account for the
  /// admixture ratio of the X chromosome
  ThetaType getTheta(const PsiType& psi) const;

  /// Return a (read-only) flattened version of the internal vectors
  const double* flat() const;

  /// Return a (read-only) flattened version of the internal vectors adjusted
  /// to account for the admixture ratio of the X chromosome
  const double* flatXChromosome(const PsiType& psi) const;

  /// Print the content of the internal vectors
  void print() const;

 private:

  /// Internal representation of the admixture proportions
  ThetaType gametes;

  /// Flat representation of the internal vectors
  ///
  /// This is updated only through calls to flat(), and should be used for
  /// read access only, as changes applied to this array are not passed back
  /// to the internal vectors.
  mutable double *vflat;

  /// Whether the dimensions internal vectors have been set
  bool isOk;

};


/** @} */


#endif /* ADMIXTUREPROPORTIONS_H */
