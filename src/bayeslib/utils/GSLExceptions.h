// *-*-C++-*-*
#ifndef __GSLEXCEPTIONS_H__
#define __GSLEXCEPTIONS_H__

//Throw exceptions instead of using GSL error handler function which
//prefers to call abort().

//Remember to put `gsl_set_error_handler(&GSLErrorHandler);' otherwise it's
//useless!

//remember also to compile with -fexceptions or they will not work!

#include <gsl/gsl_errno.h>
#include <string>
#include <exception>

using std::string;

///generic error and base struct
class GSLerror : public std::exception{   
 public:
  std::string reason, file, line;
  GSLerror(){};
  GSLerror(string r, string f, string l) : 
    reason(r), file(f), line(l) {};
  ~GSLerror() throw(){};
  const char* what()const throw(){
    return (reason + " on line " + line + " of " + file).c_str();
    //    return (file + ":" + line + " -- " + reason).c_str();
  }
};

//Structs to be thrown as exceptions

//GSL_FAILURE  = -1,

struct noConvergence ;
//GSL_CONTINUE = -2,  /* iteration has not converged */
struct badDomain ;
//GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
struct badRange ;
//GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
struct badPointer ;
//GSL_EFAULT   = 3,   /* invalid pointer */
struct badArgument ;
//GSL_EINVAL   = 4,   /* invalid argument supplied by user */
struct failure ;
//GSL_EFAILED  = 5,   /* generic failure */
struct failedFactorisation ;
//GSL_EFACTOR  = 6,   /* factorization failed */
struct failedSanity ;
//GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
struct outOfMemory ;
//GSL_ENOMEM   = 8,   /* malloc failed */
struct badFunction ;
//GSL_EBADFUNC = 9,   /* problem with user-supplied function */
struct runAway ;
//GSL_ERUNAWAY = 10,  /* iterative process is out of control */
struct maxIterations ;
//GSL_EMAXITER = 11,  /* exceeded max number of iterations */
struct divideByZero ;
//GSL_EZERODIV = 12,  /* tried to divide by zero */
struct badTolerance ;
//GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
struct aboveTolerance ;
//GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
struct underflow ;
//GSL_EUNDRFLW = 15,  /* underflow */
struct overflow ;
//GSL_EOVRFLW  = 16,  /* overflow  */
struct lossOfAccuracy ;
//GSL_ELOSS    = 17,  /* loss of accuracy */
struct roundOffError ;
//GSL_EROUND   = 18,  /* failed because of roundoff error */
struct incomformantSizes ;
//GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
struct matrixNotSquare ;
//GSL_ENOTSQR  = 20,  /* matrix not square */
struct singularityFound ;
//GSL_ESING    = 21,  /* apparent singularity detected */
struct integralOrSeriesDivergent ;
//GSL_EDIVERGE = 22,  /* integral or series is divergent */
struct badHardware ;
//GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */ 
struct notImplemented ;
//GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
struct cacheLimitExceeded ;
//GSL_ECACHE   = 25,  /* cache limit exceeded */
struct tableLimitExceeded ;
//GSL_ETABLE   = 26,  /* table limit exceeded */
struct iterationNotProgressing ;
//GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
struct jacobiansNotImprovingSolution ;  
//GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
struct cannotReachToleranceInF ;
//GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
struct cannotReachToleranceInX ;
//GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
struct cannotReachToleranceInGradient ;
//GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
struct endOfFile ;
//GSL_EOF      = 32   /* end of file */


//A few more details about the structs we're throwing as exceptions.

struct noConvergence : public GSLerror {
  noConvergence() {};
  noConvergence(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badDomain : public GSLerror {
  badDomain() {};
  badDomain(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badRange : public GSLerror {
  badRange() {};
  badRange(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badPointer : public GSLerror {
  badPointer() {};
  badPointer(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badArgument : public GSLerror {
  badArgument() {};
  badArgument(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct failure : public GSLerror {
  failure() {};
  failure(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct failedFactorisation : public GSLerror {
  failedFactorisation() {};
  failedFactorisation(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct failedSanity : public GSLerror {
  failedSanity() {};
  failedSanity(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct outOfMemory : public GSLerror {
  outOfMemory() {};
  outOfMemory(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badFunction : public GSLerror {
  badFunction() {};
  badFunction(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct runAway : public GSLerror {
  runAway() {};
  runAway(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct maxIterations : public GSLerror {
  maxIterations() {};
  maxIterations(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct divideByZero : public GSLerror {
  divideByZero() {};
  divideByZero(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badTolerance : public GSLerror {
  badTolerance() {};
  badTolerance(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct aboveTolerance : public GSLerror {
  aboveTolerance() {};
  aboveTolerance(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct underflow : public GSLerror {
  underflow() {};
  underflow(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct overflow : public GSLerror {
  overflow() {};
  overflow(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct lossOfAccuracy : public GSLerror {
  lossOfAccuracy() {};
  lossOfAccuracy(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct roundOffError : public GSLerror {
  roundOffError() {};
  roundOffError(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct incomformantSizes : public GSLerror {
  incomformantSizes() {};
  incomformantSizes(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct matrixNotSquare : public GSLerror {
  matrixNotSquare() {};
  matrixNotSquare(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct singularityFound : public GSLerror {
  singularityFound() {};
  singularityFound(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct integralOrSeriesDivergent : public GSLerror {
  integralOrSeriesDivergent() {};
  integralOrSeriesDivergent(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct badHardware : public GSLerror {
  badHardware() {};
  badHardware(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct notImplemented : public GSLerror {
  notImplemented() {};
  notImplemented(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct cacheLimitExceeded : public GSLerror {
  cacheLimitExceeded() {};
  cacheLimitExceeded(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct tableLimitExceeded : public GSLerror {
  tableLimitExceeded() {};
  tableLimitExceeded(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct iterationNotProgressing : public GSLerror {
  iterationNotProgressing() {};
  iterationNotProgressing(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct jacobiansNotImprovingSolution : public GSLerror {
  jacobiansNotImprovingSolution() {};
  jacobiansNotImprovingSolution(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};  

struct cannotReachToleranceInF : public GSLerror {
  cannotReachToleranceInF() {};
  cannotReachToleranceInF(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct cannotReachToleranceInX : public GSLerror {
  cannotReachToleranceInX() {};
  cannotReachToleranceInX(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct cannotReachToleranceInGradient : public GSLerror {
  cannotReachToleranceInGradient() {};
  cannotReachToleranceInGradient(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};

struct endOfFile : public GSLerror {
  endOfFile() {};
  endOfFile(string r, string f, string l)  : 
    GSLerror(r,f,l) {};
};


struct indexOutOfRange : public badArgument{
  size_t i,j,m,n;
  indexOutOfRange() {};
  indexOutOfRange(string r, string f, string l)  : 
    badArgument(r,f,l) {};
};


#endif //__GSLEXCEPTIONS_H__
