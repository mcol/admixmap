#include "VectorLoop.h"

VectorLoop::VectorLoop()
{
}

VectorLoop::VectorLoop(const Vector_i& NewBase):
  _Base(NewBase),
  _Count(Vector_i(NewBase.GetNumberOfElements())),
  _DecimalCount(0)
{
}

void VectorLoop::SetBase(const Vector_i& NewBase)
{
  _Base = NewBase;
  _Count.SetNumberOfElements( _Base.GetNumberOfElements() );
  _DecimalCount = 0;
}

void VectorLoop::Reset()
{
  _DecimalCount = 0;
  _Count.SetElements(0);
}

void VectorLoop::Increment(int x)
{
  for( int i=0; i<x; i++){
    int j;
    for(j=_Base.GetNumberOfElements()-1; inc_conditional(j); j--){
      _Count(j) = 0;
    }

    if( j > -1 ){
      _Count(j)++;
      _DecimalCount++;
    } else {
      _DecimalCount = 0;
    }
  }
}

// private const helper method
int VectorLoop::get_mult(int x) const
{
  int mult = 1;
  for(int i=_Base.GetNumberOfElements()-1; i>x; i--){
      mult *= _Base(i);
  }
  return mult;
}

// private const helper method
bool VectorLoop::inc_conditional(int x) const
{
  return (x>=0 && (_Count(x) == _Base(x)-1));
}

const Vector_i& VectorLoop::GetCount() const
{
  return _Count;
}

const Vector_i VectorLoop::GetCountParent(int x) const
{
  assert(x==0 || x==1);
  Vector_i parent( _Base.GetNumberOfElements() / 2 );
  for( int i = 0; i < _Base.GetNumberOfElements() / 2; i++ ){
    parent(i) = _Count( 2 * i + x );
  }
  return parent;
}

int VectorLoop::GetDecimalCount() const
{
  return _DecimalCount;
}

int VectorLoop::GetDecimalBase() const
{
  int DecimalBase = 1;
  for( int i = 0; i < _Base.GetNumberOfElements(); i++ ){
    DecimalBase *= (int)_Base( i );
  }
  return DecimalBase;
}

int VectorLoop::GetDecimal(const Vector_i& x) const
{
  int decimal = 0;
  for( int k = 0; k < _Base.GetNumberOfElements(); k++ ){
    decimal += get_mult(k) * x(k);
  }
  return decimal;
}
