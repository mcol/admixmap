//=============================================================================
//
// Copyright (C) 2007  David O'Donnell
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
/// \file HapPair.h
/// Definition of the hapPair class.
//=============================================================================

#ifndef HAPPAIR_H
#define HAPPAIR_H


#include <ostream>



/** \addtogroup base
 * @{ */


///struct to hold a pair of haplotypes, coded as integers
class hapPair
 {
 public:
   int haps[2];

   hapPair(){
     haps[0] = haps[1] = -1;
   }
   hapPair(int a, int b){
     haps[0] = a;
     haps[1] = b;
   }

   friend 
   ///for printing a happair
   std::ostream& operator<<(std::ostream& os, const hapPair &h){
     os << h.haps[0] << " " << h.haps[1];
     return os;
   }

}; 



/** @} */

#endif
