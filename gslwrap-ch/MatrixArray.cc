#include "MatrixArray.h"

MatrixArray::MatrixArray():
  SizeOfArray(1),
  TheArray(new Matrix[1])
{
  TheArray[0] = Matrix(1,1);
}
   
MatrixArray::MatrixArray( int NumberOfElements ):
  SizeOfArray(NumberOfElements),
  TheArray(0)
{
  //? should be an assertion
  if ( NumberOfElements < 1){
    SizeOfArray = 1;
  }

  TheArray = new Matrix[SizeOfArray];
  for(int i=0;i<SizeOfArray;i++){
    TheArray[i] = Matrix(1,1);
  }
}

MatrixArray::MatrixArray(const MatrixArray& M):
  SizeOfArray(M.SizeOfArray),
  TheArray(new Matrix[M.SizeOfArray])
{
  for (int i=0; i<SizeOfArray; i++){
    TheArray[i] = M.TheArray[i];
  }
}
   
MatrixArray::MatrixArray(int NumberOfElements,int NumberOfRows,int NumberOfCols):
  SizeOfArray(NumberOfElements),
  TheArray(new Matrix[NumberOfElements])
{
  for (int i = 0; i < NumberOfElements; i++ ){
    TheArray[i] = Matrix(NumberOfRows,NumberOfCols);
  }
}

MatrixArray::~MatrixArray()
{
  delete [] TheArray;
}

MatrixArray_d MatrixArray::Double()
{
   MatrixArray_d d( SizeOfArray );

   for( int i = 0; i < SizeOfArray; i++ )
      d(i) = TheArray[i].Double();
         
   return( d );
}

MatrixArray_i MatrixArray::Integer()
{
   MatrixArray_i d( SizeOfArray );

   for( int i = 0; i < SizeOfArray; i++ )
      d(i) = TheArray[i].Integer();
         
   return( d );
}

void MatrixArray::SetDimensionsOfElements( int NumberOfRows, int NumberOfCols )
{
   for (int i = 0; i < SizeOfArray; i++ )
      TheArray[i].SetNumberOfElements( NumberOfRows, NumberOfCols );
}

void MatrixArray::SetNumberOfElementsWithDimensions( int NumberOfElements, int NumberOfRows, int NumberOfCols )
{
   delete [] TheArray;

   SizeOfArray = NumberOfElements;

   TheArray = new Matrix[ NumberOfElements ];
   
   for (int i = 0; i < SizeOfArray; i++ )
      TheArray[i].SetNumberOfElements( NumberOfRows, NumberOfCols );
}

void MatrixArray::SetNumberOfElements( int NumberOfElements )
{
   SizeOfArray = NumberOfElements;
   delete [] TheArray;

   TheArray = new Matrix[ SizeOfArray ];
}

void MatrixArray::SetElements( const Matrix &TheMatrix )
{
   for (int i = 0; i < SizeOfArray; i++)
      TheArray[i] = TheMatrix;
}

void MatrixArray::SetElements( float NewValue )
{
   for (int i = 0; i < SizeOfArray; i++ )
      TheArray[i].SetElements( NewValue );
}

int MatrixArray::GetNumberOfElements() const
{
   return(SizeOfArray);
}

std::ostream& operator<< ( std::ostream& os, const MatrixArray& M)
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
                      

Matrix &MatrixArray::operator() ( int ElementNumber ) const
{
   if ( ElementNumber >= SizeOfArray )
   {
      std::cout << "WARNING: MatrixArray::operator() Element Number ";
      std::cout << ElementNumber << " > SizeOfArray " << SizeOfArray << " accessed." << std::endl;
   }

   if ( ElementNumber < 0 )
   {
      std::cout << "WARNING: MatrixArray::operator() ElementNumber ";
      std::cout << ElementNumber << " < 0 accessed." << std::endl;
   }
   
   return( TheArray[ ElementNumber ]);
}

MatrixArray &MatrixArray::operator=( const MatrixArray &m )
{
   if ( &m == this)
      return *this;

   delete [] TheArray;
   
   SizeOfArray = m.GetNumberOfElements();
   TheArray = new Matrix[ SizeOfArray ];

   for ( int i = 0; i < SizeOfArray; i++ )
      TheArray[i] = m(i);

   return *this;
}

MatrixArray &MatrixArray::operator+=( const MatrixArray &m )
{
   MatrixArray temp;
   if ( &m == this)
      return *this;

   temp = *this;
   delete [] TheArray;
   
   SizeOfArray = m.GetNumberOfElements();
   TheArray = new Matrix[ SizeOfArray ];

   for ( int i = 0; i < SizeOfArray; i++ )
      TheArray[i] = temp(i) + m(i);

   return *this;
}

MatrixArray MatrixArray::operator/( double f ) const
{
   MatrixArray temp( SizeOfArray );

   for ( int i = 0; i < SizeOfArray; i++ )
      temp(i) = TheArray[i] / f;

   return( temp );
}

void MatrixArray::RemoveElement ( int ElementNumber )
{
   if ( ElementNumber >= SizeOfArray )
   {
      std::cout << "WARNING: MatrixArray::RemoveElement Element Number";
      std::cout << ElementNumber << " > SizeOfArray" << SizeOfArray << " accessed." << std::endl;
      exit(0);
   }

   if ( ElementNumber < 0)
   {
      std::cout << "WARNING: MatrixArray::RemoveElement ElementNumber";
      std::cout << ElementNumber << " < 0 accessed." << std::endl;
      exit(0);
   }

   Matrix *TempArray;
   TempArray = new Matrix[ SizeOfArray ];
   
   for ( int i = 0; i < SizeOfArray; i++ )
      TempArray[i] = TheArray[i];
   
   SizeOfArray --;
   
   delete [] TheArray;
   TheArray = new Matrix[ SizeOfArray ];
   
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

void MatrixArray::AddElementToEnd( const Matrix &NewElement)
{
   Matrix *TempArray;
   TempArray = new Matrix[ SizeOfArray ];
   for (int i = 0; i < SizeOfArray; i++)
      TempArray[i] = TheArray[i];
   
   delete [] TheArray;
   
   TheArray = new Matrix[ (SizeOfArray + 1) ];
   for (int i = 0; i < SizeOfArray; i++)
      TheArray[i] = TempArray[i];
   
   TheArray[SizeOfArray] = NewElement;

   SizeOfArray++;
   delete [] TempArray;
}
