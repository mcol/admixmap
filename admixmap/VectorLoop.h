// *-*-C++-*-*
#ifndef VECTORLOOP_H
#define VECTORLOOP_H 1

#include "vector_i.h"

class VectorLoop
{
public:
    VectorLoop();
    VectorLoop(const Vector_i&);
    void SetBase(const Vector_i&);

    void Increment( int );
    void Reset();
    const Vector_i GetCountParent( int ) const;
    const Vector_i& GetCount() const;
    int GetDecimalCount() const;
    int GetDecimalBase() const;
    int GetDecimal(const Vector_i& x) const;

private:
    Vector_i _Base;
    Vector_i _Count;
    int _DecimalCount;
    bool inc_conditional(int) const;

    int get_mult(int) const;
};

#endif /* ! VECTORLOOP_H */
