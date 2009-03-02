// *-*-C++-*-*
#ifndef GENETICDISTANCEUNIT_H_
#define GENETICDISTANCEUNIT_H_


/** \addtogroup base
 * @{ */


// Changes to the enum here should be matched in the implementation file
//!enum for distance units
enum GeneticDistanceUnit {basepairs, kilobases, megabases, centimorgans, Morgans,
			  N_GDUS }; ///< This must remain the last tag

const char *	    gduAsString	 ( GeneticDistanceUnit u    );
GeneticDistanceUnit gduFromString( const char *	       desc );



/** @} */



#endif /*GENETICDISTANCEUNIT_H_*/
