// *-*-C++-*-*
//=============================================================================
/// \file GeneticDistanceUnit.h
/// Definition of the GeneticDistanceUnit methods.
//=============================================================================

#ifndef GENETICDISTANCEUNIT_H_
#define GENETICDISTANCEUNIT_H_


/** \addtogroup base
 * @{ */


//-----------------------------------------------------------------------------
//
// NB: Changes to the enum values here must be matched by similar changes in the
// implementation file!
//
/// Enumeration for genetic distance units.  @sa { GeneticDistance }
//
//-----------------------------------------------------------------------------

enum GeneticDistanceUnit {basepairs, kilobases, megabases, centimorgans, Morgans,
			  N_GDUS }; ///< This must remain the last tag

const char *	    gduAsString	 ( GeneticDistanceUnit u    );
GeneticDistanceUnit gduFromString( const char *	       desc );



/** @} */



#endif /*GENETICDISTANCEUNIT_H_*/
