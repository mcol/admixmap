// *-*-C++-*-*
#ifndef __GSLERRORHANDLER_H__
#define __GSLERRORHANDLER_H__

#include "bclib/bclib.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



//Throw exceptions instead of using GSL error handler function which
//prefers to call abort().

//Remember to put `gsl_set_error_handler(&GSLErrorHandler);' otherwise it's
//useless!


///Custom error handler to be used for GSL. Throw exceptions instead of
//calling abort().
void GSLErrorHandler(const char * reason, const char * file, 
			  int line, int gsl_errno); 


/** @} */

END_BCLIB_NAMESPACE

#endif //__GSLERRORHANDLER_H__
