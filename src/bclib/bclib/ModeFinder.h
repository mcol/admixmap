// *-*-C++-*-*
#ifndef MODEFINDER_H
#define MODEFINDER_H 1

#include "bclib/bclib.h"
BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///class for finding modes of functions
class ModeFinder{

 public:
  /// Newton-Raphson algorithm to find a mode of function
   static double FindModeByNewtonRaphson(double init, double* ddf, const void* const args, 
					   double (*function)(double, const void* const), 
					   double (*dfunction)(double, const void* const),
					   double (*ddfunction)(double, const void* const));
  ///uses BoundedModeSearch to find mode, first finding two points either side if unbounded
  static double SimpleModeSearch( const void* const args,
				  double (*gradient)(double, const void* const), double (*secondDeriv)(double, const void* const),
				  bool hasLowerBound, bool hasUpperBound, const double LowerBound, const double UpperBound );
     
 private:
  /// uses modified Newton-Raphson to find two points either side of the mode then runs BoundedModeSearch
  static double UnboundedModeSearch(const void* const args, double (*gradient)(double, const void* const),
				    double (*secondDeriv)(double, const void* const), bool hasLowerBound, bool hasUpperBound, 
				    const double LowerBound, const double UpperBound );
  ///searches for mode between a and b, using only gradient function
  static double BoundedModeSearch( const double a, const double b, const void* const args,
				   double (*gradient)(double, const void* const) ) ; 
};


/** @} */

END_BCLIB_NAMESPACE

#endif
