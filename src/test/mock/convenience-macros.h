#ifndef CONVENIENCEMACROS_H_
#define CONVENIENCEMACROS_H_

/** Check if pointer points to something, throw if it doesn't */
#define THROW_IF_EMPTY(VAR) if (not VAR) \
 { \
    throw CppUnit::Exception( \
    CppUnit::Message(#VAR " is not set"), \
    CPPUNIT_SOURCELINE());  }

/** Copes vectors by their pointers, used when returning by parameter */
#define RETURN_VEC_BY_PARAM(IN, OUT) \
  (OUT)->assign((IN)->begin(), (IN)->end())

#endif /*CONVENIENCEMACROS_H_*/
