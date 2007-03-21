#ifndef DEBUGMACROS_H_
#define DEBUGMACROS_H_

/** Print file, line number, variable name and variable value. */
#define SPIT(X) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X << " == " << X << endl

/** Macro for checking if a double is NaN */
#define NANHUNT(X) if (isnan(X)) { \
  SPIT(X); \
  cerr << "File: " << __FILE__ << endl; \
  cerr << "Line: " << __LINE__ << endl; \
  throw string("Nan spotted."); \
}

/** Print four variables. */
#define SPIT4(X1, X2, X3, X4) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X1 << " == " << X1 << endl \
        << #X2 << " == " << X2 << endl \
        << #X3 << " == " << X3 << endl \
        << #X4 << " == " << X4 << endl

/** Check a variable for being NaN, print three additional variables
 * that might have contributed to NaN one. */
#define NANHUNT3(X, V1, V2, V3) if (isnan(X)) { \
  SPIT4(X, V1, V2, V3); \
  throw string("Nan spotted."); \
}

/** Check if a variable is (almost) zero */
#define ZEROHUNT3(X, V1, V2, V3) \
  if (X < 1e-200) { \
    SPIT4(X, V1, V2, V3); \
    throw string("Close-zero spotted."); \
  }

#endif /*DEBUGMACROS_H_*/
#ifndef DEBUGMACROS_H_
#define DEBUGMACROS_H_

/** Print file, line number, variable name and variable value. */
#define SPIT(X) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X << " == " << X << endl

/** Macro for checking if a double is NaN */
#define NANHUNT(X) if (isnan(X)) { \
  SPIT(X); \
  cerr << "File: " << __FILE__ << endl; \
  cerr << "Line: " << __LINE__ << endl; \
  throw string("Nan spotted."); \
}

/** Print four variables. */
#define SPIT4(X1, X2, X3, X4) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X1 << " == " << X1 << endl \
        << #X2 << " == " << X2 << endl \
        << #X3 << " == " << X3 << endl \
        << #X4 << " == " << X4 << endl

/** Check a variable for being NaN, print three additional variables
 * that might have contributed to NaN one. */
#define NANHUNT3(X, V1, V2, V3) if (isnan(X)) { \
  SPIT4(X, V1, V2, V3); \
  throw string("Nan spotted."); \
}

/** Check if a variable is (almost) zero */
#define ZEROHUNT3(X, V1, V2, V3) \
  if (X < 1e-200) { \
    SPIT4(X, V1, V2, V3); \
    throw string("Close-zero spotted."); \
  }

#endif /*DEBUGMACROS_H_*/
#ifndef DEBUGMACROS_H_
#define DEBUGMACROS_H_

/** Print file, line number, variable name and variable value. */
#define SPIT(X) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X << " == " << X << endl

/** Macro for checking if a double is NaN */
#define NANHUNT(X) if (isnan(X)) { \
  SPIT(X); \
  cerr << "File: " << __FILE__ << endl; \
  cerr << "Line: " << __LINE__ << endl; \
  throw string("Nan spotted."); \
}

/** Print four variables. */
#define SPIT4(X1, X2, X3, X4) \
  cerr  << __FILE__   << ":" \
        << __LINE__   << ": " \
        << #X1 << " == " << X1 << endl \
        << #X2 << " == " << X2 << endl \
        << #X3 << " == " << X3 << endl \
        << #X4 << " == " << X4 << endl

/** Check a variable for being NaN, print three additional variables
 * that might have contributed to NaN one. */
#define NANHUNT3(X, V1, V2, V3) if (isnan(X)) { \
  SPIT4(X, V1, V2, V3); \
  throw string("Nan spotted."); \
}

/** Check if a variable is (almost) zero */
#define ZEROHUNT3(X, V1, V2, V3) \
  if (X < 1e-200) { \
    SPIT4(X, V1, V2, V3); \
    throw string("Close-zero spotted."); \
  }

#endif /*DEBUGMACROS_H_*/
