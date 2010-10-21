// *-*-C++-*-*
//=============================================================================
/// \file common.h
//=============================================================================

/// <B>Many</B> standard-library (STL, iostream, etc.) include files are
/// included from here, as well as some types defined.
///
/// @warning
/// <B>WARNING!</B> this file duplicates, in the global namespace, symbols that
/// are also defined in @a <bclib/common.h>.  Furthermore, they both use the
/// same preprocessor symbol (COMMON_H) to protect against multiple includes;
/// thus changes to either file may not propagate to any particular including
/// source file.


#ifndef COMMON_H
#define COMMON_H 1


#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>


/** \addtogroup base
 * @{ */


typedef std::vector<std::string> Vector_s; ///< std vector of strings
typedef std::vector<Vector_s>    Matrix_s; ///< std vector of std vectors

/// enum for the various regression types.
/// Mlinear = Multiple linear, Mlogistic = multiple logistic
enum RegressionType
    {
    None	,
    Linear	,
    Logistic	,
    Cox		,
    Mlinear	, ///< Multiple linear
    Mlogistic	, ///< Multiple logistic
    Multiple	,
    N_RTYPES	  ///< Keep this as the last symbol, used to count the tags
    };
extern const char * getRegressionString( RegressionType ); ///< Gives short readable description


///enum for continuous/binary datatypes
enum DataType {Continuous, Binary, CoxData};


/** @} */


#endif /* !defined COMMON_H */
