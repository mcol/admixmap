// *-*-C++-*-*
#ifndef LSHTM_MINMAX
#define LSHTM_MINMAX 1

#include <climits>
#include <cfloat>

namespace lshtm
{
  inline const signed char    lshtm_type_min(signed char    x) { return SCHAR_MIN; }
  inline const unsigned char  lshtm_type_min(unsigned char  x) { return 0;         }

  inline const signed short   lshtm_type_min(signed short   x) { return SHRT_MIN;  }
  inline const unsigned short lshtm_type_min(unsigned short x) { return 0;         }

  inline const signed int     lshtm_type_min(signed int     x) { return INT_MIN;   }
  inline const unsigned int   lshtm_type_min(unsigned int   x) { return 0;         }

  inline const signed long    lshtm_type_min(signed long    x) { return LONG_MIN;  }
  inline const unsigned long  lshtm_type_min(unsigned long  x) { return 0;         }

  inline const float  lshtm_type_min(float  x) { return FLT_MIN; }
  inline const double lshtm_type_min(double x) { return DBL_MIN; }
  inline const long double lshtm_type_min(long double x) { return LDBL_MIN; }

  inline const signed char    lshtm_type_max(signed char    x) { return SCHAR_MAX; }
  inline const unsigned char  lshtm_type_max(unsigned char  x) { return UCHAR_MAX; }

  inline const signed short   lshtm_type_max(signed short   x) { return SHRT_MAX;  }
  inline const unsigned short lshtm_type_max(unsigned short x) { return USHRT_MAX; }

  inline const signed int     lshtm_type_max(signed int     x) { return INT_MAX;   }
  inline const unsigned int   lshtm_type_max(unsigned int   x) { return UINT_MAX;  }

  inline const signed long    lshtm_type_max(signed long    x) { return LONG_MAX;  }
  inline const unsigned long  lshtm_type_max(unsigned long  x) { return ULONG_MAX; }

  inline const float  lshtm_type_max(float  x) { return FLT_MAX; }
  inline const double lshtm_type_max(double x) { return DBL_MAX; }
  inline const long double lshtm_type_max(long double x) { return LDBL_MAX; }
}
#endif /* !defined LSHTM_MINMAX */
