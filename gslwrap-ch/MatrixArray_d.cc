#include "MatrixArray_d.h"

#include "matrix_d.h"
#include "matrix.h"
#include "matrix_i.h"
#include "MatrixArray.h"
#include "MatrixArray_i.h"

MatrixArray_d::MatrixArray_d():
  SizeOfArray(1),
  TheArray(new Matrix_d[1](1,1))
{
}
   
MatrixArray_d::MatrixArray_d(int NumberOfElements):
  SizeOfArray(NumberOfElements>1 ? NumberOfElements : 1),
  TheArray(new Matrix_d[NumberOfElements>1 ? NumberOfElements : 1](1,1))
{
}

MatrixArray_d::MatrixArray_d(const MatrixArray_d& M):
  SizeOfArray(M.SizeOfArray),
  TheArray(new Matrix_d[M.SizeOfArray])
{
  for (int i=0; i<SizeOfArray; i++){
    TheArray[i] = M.TheArray[i];
  }
}
   
MatrixArray_d::MatrixArray_d(int NumberOfElements,int NumberOfRows,int NumberOfCols):
  SizeOfArray(NumberOfElements),
  TheArray(new Matrix_d[NumberOfElements](NumberOfRows,NumberOfCols))
{
//   for (int i = 0; i < NumberOfElements; i++ ){
//     TheArray[i] = Matrix_d(NumberOfRows,NumberOfCols);
//   }
}

MatrixArray_d::~MatrixArray_d()
{
  delete [] TheArray;
}

MatrixArray MatrixArray_d::Float()
{
  MatrixArray d( SizeOfArray );

  for( int i = 0; i < SizeOfArray; i++ )
    d(i) = TheArray[i].Float();
         
  return( d );
}

MatrixArray_i MatrixArray_d::Integer()
{
  MatrixArray_i d( SizeOfArray );

  for( int i = 0; i < SizeOfArray; i++ )
    d(i) = TheArray[i].Integer();
         
  return( d );
}

void MatrixArray_d::SetDimensionsOfElements( int NumberOfRows, int NumberOfCols )
{
  for (int i = 0; i < SizeOfArray; i++ )
    TheArray[i].SetNumberOfElements( NumberOfRows, NumberOfCols );
}

void MatrixArray_d::SetNumberOfElementsWithDimensions( int NumberOfElements, int NumberOfRows, int NumberOfCols )
{
  delete [] TheArray;

  SizeOfArray = NumberOfElements;

  TheArray = new Matrix_d[ NumberOfElements ];
   
  for (int i = 0; i < SizeOfArray; i++ )
    TheArray[i].SetNumberOfElements( NumberOfRows, NumberOfCols );
}

void MatrixArray_d::SetNumberOfElements( int NumberOfElements )
{
  SizeOfArray = NumberOfElements;
  delete [] TheArray;

  TheArray = new Matrix_d[ SizeOfArray ];
}

void MatrixArray_d::SetElements( const Matrix_d &TheMatrix )
{
  for (int i = 0; i < SizeOfArray; i++)
    TheArray[i] = TheMatrix;
}

void MatrixArray_d::SetElements( float NewValue )
{
  for (int i = 0; i < SizeOfArray; i++ )
    TheArray[i].SetElements( NewValue );
}

int MatrixArray_d::GetNumberOfElements() const
{
  return(SizeOfArray);
}

std::ostream& operator<< ( std::ostream& os, const MatrixArray_d& M)
{
  int i;

  os.setf( std::ios::fixed );

  for ( i = 0; i < M.GetNumberOfElements(); i++ )
    {
      os << "[" << i << "]" << std::endl;
      os << M(i) << std::endl;
    }

  return os;
}
                      

Matrix_d &MatrixArray_d::operator() ( int ElementNumber ) const
{
   assert( ElementNumber < SizeOfArray );
   assert( ElementNumber >= 0 );
   
  return( TheArray[ ElementNumber ]);
}

MatrixArray_d &MatrixArray_d::operator=( const MatrixArray_d &m )
{
  if ( &m == this)
    return *this;

  delete [] TheArray;
   
  SizeOfArray = m.GetNumberOfElements();
  TheArray = new Matrix_d[ SizeOfArray ];

  for ( int i = 0; i < SizeOfArray; i++ )
    TheArray[i] = m(i);

  return *this;
}

MatrixArray_d &MatrixArray_d::operator+=( const MatrixArray_d &m )
{
  MatrixArray_d temp;
  if ( &m == this)
    return *this;

  temp = *this;
  delete [] TheArray;
   
  SizeOfArray = m.GetNumberOfElements();
  TheArray = new Matrix_d[ SizeOfArray ];

  for ( int i = 0; i < SizeOfArray; i++ )
    TheArray[i] = temp(i) + m(i);

  return *this;
}

MatrixArray_d MatrixArray_d::operator/( double f ) const
{
  MatrixArray_d temp( SizeOfArray );

  for ( int i = 0; i < SizeOfArray; i++ )
    temp(i) = TheArray[i] / f;

  return( temp );
}

void MatrixArray_d::RemoveElement ( int ElementNumber )
{
  if ( ElementNumber >= SizeOfArray )
    {
      std::cout << "WARNING: MatrixArray_d::RemoveElement Element Number ";
      std::cout << ElementNumber << " > SizeOfArray " << SizeOfArray << " accessed." << std::endl;
      exit(0);
    }

  if ( ElementNumber < 0)
    {
      std::cout << "WARNING: MatrixArray_d::RemoveElement ElementNumber";
      std::cout << ElementNumber << " < 0 accessed." << std::endl;
      exit(0);
    }

  Matrix_d *TempArray;
  TempArray = new Matrix_d[ SizeOfArray ];
   
  for ( int i = 0; i < SizeOfArray; i++ )
    TempArray[i] = TheArray[i];
   
  SizeOfArray --;
   
  delete [] TheArray;
  TheArray = new Matrix_d[ SizeOfArray ];
   
  int counter = -1;
  for ( int i = 0; i < (SizeOfArray + 1); i++ )
    {
      if (i != ElementNumber)
	{
	  counter++;
	  TheArray[counter] = TempArray[i];
	}
    }
  delete [] TempArray;
}

void MatrixArray_d::AddElementToEnd( const Matrix_d &NewElement)
{
  Matrix_d *TempArray;
  TempArray = new Matrix_d[ SizeOfArray ];
  for (int i = 0; i < SizeOfArray; i++)
    TempArray[i] = TheArray[i];
   
  delete [] TheArray;
   
  TheArray = new Matrix_d[ (SizeOfArray + 1) ];
  for (int i = 0; i < SizeOfArray; i++)
    TheArray[i] = TempArray[i];
   
  TheArray[SizeOfArray] = NewElement;

  SizeOfArray++;
  delete [] TempArray;
}

