//=============================================================================
//
// Copyright (C) 2009  David D. Favro  amm-header@meta-dynamic.com
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 3 as published by the Free
// Software Foundation.
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
/// \file NRand.h
/// Interface of the NRand class.
//=============================================================================


#ifndef __bclib_bclib_NRand_h
#define __bclib_bclib_NRand_h



#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



namespace genepi { // ----
//namespace bclib { // ----

/** \addtogroup bclib
 * @{ */



//-----------------------------------------------------------------------------
/// Random-number generator: essentially a C++ wrapper around the GSL C routines
/// (see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html),
/// although it could be reimplemented with a different backend.  While intended
/// to replace Rand (which, amongst other problems, is not thread-safe), it
/// offers some GSL functions that Rand lacks, and lacks some of the
/// functionality that Rand provides.  The meaning of const-ness here is rather
/// arbitrary: we assume changing the internal state of the PRNG
/// (i.e. generating or seeding) to require a non-const reference.
//-----------------------------------------------------------------------------

class NRand
    {
    private:
	gsl_rng * const gen;

    public:
	NRand( const gsl_rng_type * type = gsl_rng_taus, unsigned long init_seed = gsl_rng_default_seed );
	NRand( unsigned long init_seed ); ///< Allocates a Tausworthe generator

	~NRand();

	void seed( unsigned long new_seed );

	/// The raw PRNG.
	unsigned long get()	  { return gsl_rng_get( gen ); }
	unsigned long min() const { return gsl_rng_min( gen ); }
	unsigned long max() const { return gsl_rng_max( gen ); }

	double standard_uniform	   () { return gsl_rng_uniform	  ( gen ); } ///< [0,1)
	double standard_uniform_pos() { return gsl_rng_uniform_pos( gen ); } ///< (0,1)

	unsigned long uniform_int( unsigned long n ) { return gsl_rng_uniform_int( gen, n ); } ///< [0,n-1]

	double uniform( double min, double max );
	double normal( double mean, double sigma ) { return mean + gsl_ran_gaussian( gen, sigma ); }
    };


/** @} */

//} // ---- end namespace bclib
} // ---- end namespace genepi



#endif // ! __bclib_bclib_NRand_h
