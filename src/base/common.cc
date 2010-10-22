// *-*-C++-*-*
//=============================================================================
/// \file common.cc
/// Implementation of the getRegressionString(RegressionType) function.
//=============================================================================

#include "common.h"

#include <stdexcept>


static const char * const RegressionString [] =
    { "None", "Linear", "Logistic", "Cox", "Mlinear", "Mlogistic", "Multiple" };

const char * getRegressionString( RegressionType t )
    {
    if ( (sizeof(RegressionString)/sizeof(*RegressionString)) != N_RTYPES )
	throw std::range_error( "failed N_TYPES assertion" );

    // STL is so bad: we want: std::string("Regression-type ") + int(t) " is out-of-bounds" );
    if ( (int(t) < 0) || (int(t) >= N_RTYPES) )
	throw std::invalid_argument( "Regression-type is out-of-bounds" );

    return RegressionString[ t ];
    }
