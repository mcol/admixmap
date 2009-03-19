//=============================================================================
//
// Copyright (C) 2009  David D. Favro  gpl-copyright@meta-dynamic.com
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
/// \file AncestryVector.cc
/// Implementation of the AncestryVector class
//=============================================================================

#include "AncestryVector.h"


#if AV_OSTREAM
    #include <iostream>
#endif



namespace genepi { // ----



//-----------------------------------------------------------------------------
// to_ulong()
//-----------------------------------------------------------------------------

unsigned long AncestryVector::to_ulong() const
    {
    unsigned long rv = 0;

    if ( K == 2 )	// Optimized version for K==2:
	for ( size_t idx = size() ; idx-- != 0 ; )
	    rv = (rv << 1) + at(idx);

    else if ( K == 3 )	// Optimized version for K==3:
	for ( size_t idx = size() ; idx-- != 0 ; )
	    rv = (rv << 1) + rv + at(idx);

    else
	for ( size_t idx = size() ; idx-- != 0 ; )
	    rv = (rv * K) + at(idx);


    // We need to cache the limit here:
    //gp_assert_lt( rv, K^F );


    return rv;
    }



//-----------------------------------------------------------------------------
// Print an AncestryVector to an ostream
//-----------------------------------------------------------------------------

#if AV_OSTREAM

    static AVOutputStyle style = AV_NUMERIC;

    void setAVOutputStyle( AVOutputStyle nv )
	{
	style = nv;
	}

    inline static std::ostream & output( std::ostream & os, PopIdx pop )
	{
	#if 0
	    if ( style == AV_NUMERIC )
		return os << popVect[pop];
	    else
	#endif
		return os << size_t(pop);
	}

    std::ostream & operator<<( std::ostream & os, const AncestryVector & av )
	{
	gp_assert( av.size() != 0 );

	os << "AV(";
	const size_t limit = av.size() - 1;
	for ( size_t idx = 0 ; idx < limit ; ++idx )
	    output( os, av.at(idx) ) << ',';
	os << av.at(limit) << ')';

	return os;
	}

#endif // AV_OSTREAM



} // ---- end namespace genepi
