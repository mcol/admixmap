#include "bclib/pvector.h"
#include <cmath>//for fabs
#include <algorithm>
#include <numeric>
//#include <string>//for exceptions
#include "IsNegative.h" //for verify
#include <iterator>


#include <gsl/gsl_sf_log.h> // For softmax()
#include <gsl/gsl_sys.h>    // For softmax() [gsl_finite()]
#include <gsl/gsl_errno.h>  // For softmax()
#include <sstream>	    // For softmax()
#include <algorithm>	    // For softmax() [max()]


#define PVECTOR_PRECISION 1e-10
#define PVECTOR_SUM 1.0

BEGIN_BCLIB_NAMESPACE

template <class T>
bool pvector<T>::is_normalized(){
  // TODO: Accumulate only when necessary, i.e.
  // when vector is changed.
  sum = accumulate(this->begin(), this->end(), 0.0);
  return (fabs(sum - PVECTOR_SUM) < PVECTOR_PRECISION);
}

template <class T>
void pvector<T>::normalize(){
  if (!this->is_normalized()) { // Updates sum
    (*this) *= 1.0 / sum; 
  }
}

template <class T>
bool pvector<T>::verify(){
  if (!this->is_normalized()) {
    throw ("pvector<T>::verify(): doesn't sum to one.");
  }
  if((unsigned)count_if(this->begin(), this->end(), IsNegative<T>) > 1) {
    throw ("pvector<T>::verify(): Some elements are less than zero.");
  }
  return true;
}

template <class T>
void pvector<T>::snapToZero(){
 for(typename pvector<T>::iterator i = this->begin(); i != this->end(); ++i){
  if(*i < threshold) *i = 0;
 } 
 //snapToZero(threshold); // DDF: WTF: This appears to cause infinite recursion, so I commented it out.
 this->normalize();
}



//=============================================================================
// Copy-result-out (not in-place) transformations
//=============================================================================

template < typename T > template< typename DestIter > void pvector<T>::inv_softmax( DestIter dest ) const
    {
    gsl_sf_result result;

    // Disable default gsl error handler:
    gsl_error_handler_t * const old_handler = gsl_set_error_handler_off();

    T logz = 0.0;

    DestIter destPtr = dest;
    for ( const_iterator it = begin(); it != end(); ++it )
	{
	const int status = gsl_sf_log_e( *it, &result );
	if ( status != 0 )
	    {
	    std::stringstream err;
	    err << "error in inv_softmax: " << gsl_strerror(status) << ": " << *it;
	    throw err.str();
	    }
	*destPtr++ = result.val;
	logz -= result.val;
	}

    logz /= T( size() );

    for ( const_iterator it = begin(); it != end(); ++it )
	{
	const T newVal = *dest + logz;
	*dest++ = newVal;
	if ( ! gsl_finite( newVal ) )
	    throw std::string( "error in inv_softmax" );
	}

    gsl_set_error_handler( old_handler ); //restore gsl error handler
    }



template < typename T > template< typename DestIter, typename QualFunctor >
		void pvector<T>::inv_softmax( DestIter dest, QualFunctor qualifies, const T & defVal ) const
    {

    gsl_sf_result result;

    // Disable default gsl error handler:
    gsl_error_handler_t * const old_handler = gsl_set_error_handler_off();

    bool qualArray[ size() ];

    T logz = 0.0;

    DestIter destPtr = dest;
    bool *   qualPtr = qualArray;
    for ( const_iterator it = begin(); it != end(); ++it, ++destPtr )
	{
	const T	   val = *it;
	const bool qual = qualifies( val );
	*qualPtr++ = qual;
	if ( qual )
	    {
	    const int status = gsl_sf_log_e( *it, &result );
	    if ( status != 0 )
		{
		std::stringstream err;
		err << "error in inv_softmax: " << gsl_strerror(status) << ": " << *it;
		throw err.str();
		}
	    *destPtr = result.val;
	    logz -= result.val;
	    }
	else
	    *destPtr = defVal;
	}

    // Perhaps this should divide by the number of _qualifying_ elements?
    logz /= T( size() );

    qualPtr = qualArray;
    for ( const_iterator it = begin(); it != end(); ++it, ++dest )
	if ( *qualPtr++ )
	    {
	    const T newVal = *dest + logz;
	    *dest = newVal;
	    if ( ! gsl_finite( newVal ) )
		throw std::string( "error in inv_softmax" );
	    }

    gsl_set_error_handler( old_handler ); //restore gsl error handler

    }



template < typename T > template< typename DestIter, typename QualFunctor >
		void pvector<T>::softmax( DestIter dest, QualFunctor qualifies, const T & defVal ) const
    {

    // Check for null input array, which will cause a singularity:
    if ( size() == 0 )
	throw std::runtime_error( "Null input to softmax()" );


    bool qualArray[ size() ];

    bool * qualPtr = qualArray;

    T amax = 0.0;
    for ( const_iterator it = begin(); it != end(); ++it )
	{
	const T val = *it;
	const bool qual = qualifies( val );
	*qualPtr++ = qual;
	if ( qual )
	    amax = std::max( amax, *it );
	}


    T total = 0.0;

    DestIter destPtr = dest;
    qualPtr = qualArray;
    for ( const_iterator it = begin(); it != end(); ++it, ++destPtr )
	if ( *qualPtr++ )
	    {
	    const T val = exp( *it - amax );
	    total += val;
	    *destPtr = val;
	    }
	else
	    *destPtr = defVal;


    qualPtr = qualArray;
    for ( const_iterator it = begin(); it != end(); ++it, ++dest )
	if ( *qualPtr++ )
	    *dest /= total;

    }



//=============================================================================
// Providing explicit template instantiation to avoid linkage errors.
//=============================================================================

template bool pvector<double>::is_normalized();
template void pvector<double>::normalize();
template bool pvector<double>::verify();
template void pvector<double>::snapToZero();

template bool pvector<float>::is_normalized();
template void pvector<float>::normalize();
template bool pvector<float>::verify();
template void pvector<float>::snapToZero();


template void pvector<double>::inv_softmax<double*>( double* ) const;
template void pvector<double>::inv_softmax<genepi::cvector<double>::iterator>( pvector<double>::iterator ) const;

typedef bool (* DoubleTester)(double);
template void pvector<double>::inv_softmax<double*,DoubleTester>( double*, DoubleTester, const double& ) const;
template void pvector<double>::inv_softmax<genepi::cvector<double>::iterator,DoubleTester>( pvector<double>::iterator, DoubleTester, const double& ) const;

template void pvector<double>::inv_softmax_gt0<double*>( double*, const double& ) const;
template void pvector<double>::inv_softmax_gt0<genepi::cvector<double>::iterator>( pvector<double>::iterator, const double& ) const;

template void pvector<double>::softmax<pvector<double>::iterator,DoubleTester>( genepi::cvector<double>::iterator, DoubleTester, const double& ) const;
template void pvector<double>::softmax<double*,DoubleTester>( double*, DoubleTester, const double & ) const;



END_BCLIB_NAMESPACE
