#include "MatrixArray_i.h"

#include "matrix_i.h"
#include "matrix_d.h"
#include "matrix.h"
#include "MatrixArray_d.h"
#include "MatrixArray.h"

MatrixArray_i::MatrixArray_i():
  SizeOfArray(1),
  TheArray(new Matrix_i[1])
{
  Matrix_i nullMatrix(1,1);
  TheArray[0] = nullMatrix;
}
   
MatrixArray_i::MatrixArray_i(int NumberOfElements):
  SizeOfArray(NumberOfElements),
  TheArray(0)
{
  //? should be assertion
  if (NumberOfElements < 1){
    SizeOfArray = 1;
  }

  Matrix_i nullMatrix(1,1);
  TheArray = new Matrix_i[ SizeOfArray ];
  for(int i=0;i<SizeOfArray;i++){
    TheArray[i] = nullMatrix;
  }
}

MatrixArray_i::MatrixArray_i(const MatrixArray_i& M):
  SizeOfArray(M.SizeOfArray),
  TheArray(new Matrix_i[M.SizeOfArray])
{
  for (int i = 0; i<SizeOfArray; i++){
    TheArray[i] = M.TheArray[i];
  }
}
   
MatrixArray_i::MatrixArray_i(int NumberOfElements,int NumberOfRows,int NumberOfCols):
  SizeOfArray(NumberOfElements),
  TheArray(new Matrix_i[NumberOfElements])
{
  for (int i=0; i<NumberOfElements; i++){
    TheArray[i] = Matrix_i(NumberOfRows,NumberOfCols);
  }
}

MatrixArray_i::~MatrixArray_i()
{
  delete [] TheArray;
}

MatrixArray MatrixArray_i::Float()
{
  MatrixArray d( SizeOfArray );

  for( int i = 0; i < SizeOfArray; i++ )
    d(i) = TheArray[i].Float();
         
  return( d );
}

MatrixArray_d MatrixArray_i::Double()
{
  MatrixArray_d d( SizeOfArray );

  for( int i = 0; i < SizeOfArray; i++ )
    d(i) = TheArray[i].Double();
         
  return( d );
}

void MatrixArray_i::SetDimensionsOfElements( int NumberOfRows, int NumberOfCols )
{
  for (int i = 0; i < SizeOfArray; i++ )
    TheArray[i].SetNumberOfElements( NumberOfRows, NumberOfCols );
}

void MatrixArray_i::SetNumberOfElementsWithDimensions( int NumberOfElements, int NumberOfRows, int NumberOfCols )
{
  delete [] TheArray;

  SizeOfArray = NumberOfElements;

  TheArray = new Matrix_i[ NumberOfElements ];
   
  for (int i = 0; i < SizeOfArray; i++ )
    TheArray[i].SetNumberOfElements( NumberOfRows, NumberOfCols );
}

void MatrixArray_i::SetNumberOfElements( int NumberOfElements )
{
  SizeOfArray = NumberOfElements;
  delete [] TheArray;

  TheArray = new Matrix_i[ SizeOfArray ];
}

void MatrixArray_i::SetElements( const Matrix_i &TheMatrix )
{
  for (int i = 0; i < SizeOfArray; i++)
    TheArray[i] = TheMatrix;
}

void MatrixArray_i::SetElements( int NewValue )
{
  for (int i = 0; i < SizeOfArray; i++ )
    TheArray[i].SetElements( NewValue );
}

int MatrixArray_i::GetNumberOfElements() const
{
  return(SizeOfArray);
}

std::ostream& operator<< ( std::ostream& os, const MatrixArray_i& M)
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
                      

Matrix_i &MatrixArray_i::operator() ( int ElementNumber ) const
{
   assert( ElementNumber < SizeOfArray );
   assert( ElementNumber >= 0 );
//   if ( ElementNumber >= SizeOfArray )
//     {
//       std::cout << "WARNING: MatrixArray_i::operator() Element Number";
//       std::cout << ElementNumber << " > SizeOfArray" << SizeOfArray << " accessed." << std::endl;
//       exit(0);
//     }

//   if ( ElementNumber < 0)
//     {
//       std::cout << "WARNING: MatrixArray_i::operator() ElementNumber";
//       std::cout << ElementNumber << " < 0 accessed." << std::endl;
//       exit(0);
//     }
   
  return( TheArray[ ElementNumber ]);
}

MatrixArray_i &MatrixArray_i::operator=( const MatrixArray_i &m )
{
  if ( &m == this)
    return *this;

  delete [] TheArray;
   
  SizeOfArray = m.GetNumberOfElements();
  TheArray = new Matrix_i[ SizeOfArray ];

  for ( int i = 0; i < SizeOfArray; i++ )
    TheArray[i] = m(i);

  return *this;
}

MatrixArray_i &MatrixArray_i::operator+=( const MatrixArray_i &m )
{
  MatrixArray_i temp;
  if ( &m == this)
    return *this;

  temp = *this;
  delete [] TheArray;
   
  SizeOfArray = m.GetNumberOfElements();
  TheArray = new Matrix_i[ SizeOfArray ];

  for ( int i = 0; i < SizeOfArray; i++ )
    TheArray[i] = temp(i) + m(i);

  return *this;
}

void MatrixArray_i::RemoveElement ( int ElementNumber )
{
  if ( ElementNumber >= SizeOfArray )
    {
      std::cout << "WARNING: MatrixArray_i::RemoveElement Element Number";
      std::cout << ElementNumber << " > SizeOfArray" << SizeOfArray << " accessed." << std::endl;
      exit(0);
    }

  if ( ElementNumber < 0)
    {
      std::cout << "WARNING: MatrixArray_i::RemoveElement ElementNumber";
      std::cout << ElementNumber << " < 0 accessed." << std::endl;
      exit(0);
    }

  Matrix_i *TempArray;
  TempArray = new Matrix_i[ SizeOfArray ];
   
  for ( int i = 0; i < SizeOfArray; i++ )
    TempArray[i] = TheArray[i];
   
  SizeOfArray --;
   
  delete [] TheArray;
  TheArray = new Matrix_i[ SizeOfArray ];
   
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

void MatrixArray_i::AddElementToEnd( const Matrix_i &NewElement)
{
  Matrix_i *TempArray;
  TempArray = new Matrix_i[ SizeOfArray ];
  for (int i = 0; i < SizeOfArray; i++)
    TempArray[i] = TheArray[i];
   
  delete [] TheArray;
   
  TheArray = new Matrix_i[ (SizeOfArray + 1) ];
  for (int i = 0; i < SizeOfArray; i++)
    TheArray[i] = TempArray[i];
   
  TheArray[SizeOfArray] = NewElement;

  SizeOfArray++;
  delete [] TempArray;
}

