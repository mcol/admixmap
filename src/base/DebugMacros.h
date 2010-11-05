//=============================================================================
//
// Copyright (C) 2007  Maciej Bliziński
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file DebugMacros.h
/// Definition of some debugging macros.
//=============================================================================

#ifndef DEBUGMACROS_H_
#define DEBUGMACROS_H_

/** Print file, line number, variable name and variable value. */
#define SPIT(X) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X << " == " << X << endl

/** Macro for checking if a double is NaN */
#define NANHUNT(X) if (isnan(X)) { \
  SPIT(X); \
  cerr << "File: " << __FILE__ << endl; \
  cerr << "Line: " << __LINE__ << endl; \
  throw string("Nan spotted."); \
}

/** Print four variables. */
#define SPIT4(X1, X2, X3, X4) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X1 << " == " << X1 << endl \
        << #X2 << " == " << X2 << endl \
        << #X3 << " == " << X3 << endl \
        << #X4 << " == " << X4 << endl

/** Check a variable for being NaN, print three additional variables
 * that might have contributed to NaN one. */
#define NANHUNT3(X, V1, V2, V3) if (isnan(X)) { \
  SPIT4(X, V1, V2, V3); \
  throw string("Nan spotted."); \
}

/** Check if a variable is (almost) zero */
#define ZEROHUNT3(X, V1, V2, V3) \
  if (X < 1e-200) { \
    SPIT4(X, V1, V2, V3); \
    throw string("Close-zero spotted."); \
  }
  
#define CONDITIONAL_SPIT3(CONDITION, V0, V1, V2) \
  if (CONDITION) { \
    SPIT(V0); \
    SPIT(V1); \
    SPIT(V2); \
  }


#endif /*DEBUGMACROS_H_*/
