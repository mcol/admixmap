#include "vector_d.h"

#include "rand.h"
#include "matrix_d.h"
#include "vector_i.h"
#include "matrix_i.h"
#include "vector.h"
#include "matrix.h"

#define BOUNDS_CHECK 1

using namespace std;

Vector_d::Vector_d():
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
}

Vector_d::Vector_d(const Vector_d& V):
  NumberOfElements(V.NumberOfElements),
  NumberOfMissingValues(V.NumberOfMissingValues),
  MissingValues(0),
  vector(gsl_vector_calloc(V.NumberOfElements))
{

  // Copy vector elements
  gsl_vector_memcpy( vector, V.vector );

  // Copy missing values
  if (NumberOfMissingValues){
    int i;
    MissingValues = new int[NumberOfMissingValues];
    for (i=0; i<NumberOfMissingValues; i++){
      MissingValues[i] = V.MissingValues[i];
    }
  }
}

Vector_d::Vector_d(int NewNumberOfElements):
  NumberOfElements(NewNumberOfElements),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  //? should be an assertion
  if ( NewNumberOfElements <= 0 ){
    NumberOfElements = 1;
  }

  vector = gsl_vector_calloc( NumberOfElements );
}

Vector_d::~Vector_d()
{
  if (vector != 0){
    gsl_vector_free(vector);
  }
  if (NumberOfMissingValues){
    delete [] MissingValues;
  }
}

Vector Vector_d::Float()
{
  Vector d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (float)(gsl_vector_get( vector, i ) );
   
  return( d );
}
                      
double &Vector_d::operator()( int element )
{
#ifdef BOUNDS_CHECK
   assert( element < NumberOfElements );
   assert( element > -1 );
//   if ( element > NumberOfElements - 1 )
//     {
//       cout << "index > NumberOfElements - 1 in Vector_d::().\n";
//       exit(0);
//     }
//   if ( element < 0 )
//     {
//       cout << "index < 0 in Vector_d::().\n";
//       exit(0);
//     }
#endif
  return( vector->data[element*vector->stride] );
}

double& Vector_d::operator()(int i) const{
#ifdef BOUNDS_CHECK   
   assert( i < NumberOfElements );
   assert( i > -1 );
//   if ( i > NumberOfElements - 1 ){
//     cout << "index = " << i << " > NumberOfElements - 1 in Vector_d::().\n";
//     exit(0);
//   }
//   if ( i < 0 ){
//     cout << "index = " << i << " < 0 in Vector_i::().\n";
//     exit(0);
//   }
#endif
  return( vector->data[i*vector->stride] );  
}

int Vector_d::GetNumberOfElements() const
{
  return( NumberOfElements );
}

void Vector_d::SetNumberOfElements( int NewValue )
{
  if ( NumberOfElements == NewValue )
    {
      return;
    }

  if ( vector != 0 ){
    gsl_vector_free( vector );
  }

  if ( NumberOfMissingValues != 0 )
    {
      delete [] MissingValues;
      NumberOfMissingValues = 0;
      MissingValues = 0;
    }

  if ( NewValue <= 0 )
    {
      NewValue = 1;
    }

  NumberOfElements = NewValue;
  vector = gsl_vector_calloc( NumberOfElements );
}

double Vector_d::GetElement( int index ) const
{
#ifdef BOUNDS_CHECK   
   assert( index < NumberOfElements );
   assert( index > -1 );
//   if ( index > NumberOfElements - 1 )
//     {
//       cout << "index > NumberOfElements - 1 in Vector_d::GetElement.\n";
//       exit(0);
//     }
//   if ( index < 0 )
//     {
//       cout << "index < 0 in Vector_d::GetElement.\n";
//       exit(0);
//     }
#endif
  return( (double)gsl_vector_get( vector, index ) );
}

ostream& operator<< ( ostream& os, const Vector_d& V )
{
  int i;

  for ( i = 0; i < V.GetNumberOfElements() - 1; i++ )
    {
      if ( V.IsMissingElement( i ) )
	os << "# ";
      else
	os << V.GetElement( i ) << " ";
    }
  if ( V.IsMissingElement( i ) )
    os << "#" << endl;
  else
    os << V.GetElement( i );

  return os;
}	

istream& operator>> ( istream& is, Vector_d& V )
{
  int i;
  int NumCols;
  int *TempMissingValues;
  char *CurrentLine;
  char delimiters[] = " \t";
  char *p;
  char PreviousChar;
  char *token;
  char *TailPointer;

  CurrentLine = new char[LINELENGTH];
   
  // Save initial stream position and scan through file counting rows
  // and columns

  is.getline( CurrentLine, LINELENGTH );
  p = CurrentLine;
  PreviousChar = *p;
  NumCols = 0;

  while ( *p )
    {
      // Previous char not a delimiter and current character is a
      // delimiter
       
      if ( strchr( delimiters, PreviousChar ) == NULL && strchr( delimiters, *p ) != NULL )
	NumCols++;
      PreviousChar = *p;
      p++;
    }

  // If final character before EOL was not a delimiter add a final column

  if ( strchr( delimiters, PreviousChar ) == NULL )
    NumCols++;

  // Last line of file contains no columns, so suck in chars until eof
   
  if ( NumCols == 0 )
    {
      cout << "ERROR: Line contains no columns." << endl;
      exit( 1 );
    }

  // Allocate vector of the correct dimensions
   
  V.SetNumberOfElements( NumCols );

  // Worst case is that all values are missing
   
  TempMissingValues = new int[NumCols];

  // Break up tokens using delimiters and strtok
  // Missing values are coded as '#'
   
  i = 0;
  V.NumberOfMissingValues = 0;
  while ( ( token = strtok ( ( i == 0 ) ? CurrentLine : NULL, delimiters ) ) )
    {
      if ( !strcmp( token, "#" ) ){
	V( i ) = 0.0;
	TempMissingValues[V.NumberOfMissingValues++] = i;
      } else {
	V( i ) = strtod( token, &TailPointer );
      }
      i++;
    }

  // Copy missing value indices from TempMissingValues to V.MissingValues
   
  if ( V.NumberOfMissingValues )
    V.MissingValues = new int[V.NumberOfMissingValues];
  else
    V.MissingValues = 0;
  for ( i = 0; i < V.NumberOfMissingValues; i++ )
    V.MissingValues[i] = TempMissingValues[i];

  // Peek one character ahead to test for eof, put char back if we're not at eof
   
  is.get( PreviousChar );
  if ( is )
    is.putback( PreviousChar );
   
  delete [] CurrentLine;
  delete [] TempMissingValues;
   
  return is;
}	

Vector_d Vector_d::operator+( const Vector_d &v ) const
{
  int i;
  int offset;
  Vector_d temp( NumberOfElements );

  assert( v.NumberOfElements == NumberOfElements );
//     {
//       cout << "ERROR: Attempt to add vectors(_d) of different length." << endl;
//       exit( 1 );
//     }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + v.vector->data[offset];
    }

  return( temp );
}

Vector_d Vector_d::operator+( const Vector_i &v ) const
{
  int i;
  int offset;
  Vector_d temp( NumberOfElements );

  assert( v.GetNumberOfElements() == NumberOfElements );
//     {
//       cout << "ERROR: Attempt to add vectors(_d) of different length." << endl;
//       exit( 1 );
//     }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + (int)v(i);
    }

  return( temp );
}

Vector_d &Vector_d::operator+=( const Vector_d &v )
{
  int i;
  int offset;

  assert( v.NumberOfElements == NumberOfElements );
//     {
//       cout << "ERROR: Attempt to add vectors(_d) of different length." << endl;
//       exit( 1 );
//     }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] += v.vector->data[offset];
    }

  return( *this );
}

Vector_d Vector_d::operator-( const Vector_d &v ) const
{
  int i;
  int offset;
  Vector_d temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to subtract vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] - v.vector->data[offset];
    }

  return( temp );
}

Vector_d Vector_d::operator*( double f ) const
{
  int i;
  int offset;
  Vector_d temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] * f;
    }

  return( temp );
}

Vector_d Vector_d::operator/( const Vector_d &v ) const
{
  int i;
  int offset;
  Vector_d temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to divide vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      if ( v.vector->data[offset] == 0.0 )
	{
	  cout << "ERROR: Division by zero in element " << i << " of denominator vector." << endl;
	  exit( 1 );
	}
      temp.vector->data[offset] = vector->data[offset] / v.vector->data[offset];
    }

  return( temp );
}

Vector_d Vector_d::operator/( double f ) const
{
  int i;
  int offset;
  Vector_d temp( NumberOfElements );

  assert( f!= 0.0 );
//   if ( f == 0.0 )
//     {
//       cout << "ERROR: Attempted to divide vector_d by zero." << endl;
//       exit( 1 );
//     }
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] / f;
    }

  return( temp );
}

Vector_d &Vector_d::operator/=( double f )
{
  int i;
  assert( f != 0.0 );   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] /= f;
    }

  return( *this );
}

Vector_d &Vector_d::operator=( const Vector_d &v )
{
  int i;

  // Check for self-assignment
   
  if ( &v == this )
    return *this;

  // De-allocate old vector, allocate new one
   
  if ( vector != 0 )
    gsl_vector_free( vector );
  NumberOfElements = v.GetNumberOfElements();
  vector = gsl_vector_calloc( NumberOfElements );

  // Same for missing values
   
  if ( NumberOfMissingValues )
    delete [] MissingValues;
  NumberOfMissingValues = v.GetNumberOfMissingValues();
  if ( NumberOfMissingValues )
    MissingValues = new int[NumberOfMissingValues];
  else
    MissingValues = 0;
  for ( i = 0; i < NumberOfMissingValues; i++ )
    MissingValues[i] = v.GetMissingValueIndex( i );

  // Copy contents of v into *this
  gsl_vector_memcpy( vector, v.vector );

  return *this;
}

int Vector_d::GetNumberOfMissingValues() const
{
  return( NumberOfMissingValues );
}

int Vector_d::GetMissingValueIndex( int MissingValueNumber ) const
{
#ifdef BOUNDS_CHECK   
  if ( MissingValueNumber > NumberOfMissingValues - 1 )
    MissingValueNumber = NumberOfMissingValues - 1;
  if ( MissingValueNumber < 0 )
    MissingValueNumber = 0;
#endif

  return( MissingValues[MissingValueNumber] );
}

int Vector_d::IsMissingElement( int index ) const
{
#ifdef BOUNDS_CHECK   
  if ( index > NumberOfElements - 1 )
    {
      index = NumberOfElements - 1;
    }
  if ( index < 0 )
    {
      index = 0;
    }
#endif
  int i;

  for ( i = 0; i < NumberOfMissingValues; i++ )
    if ( index == MissingValues[i] )
      return( 1 );
  return( 0 );
}

Matrix_d Vector_d::ColumnMatrix() const
{
  int i;
  Matrix_d temp( NumberOfElements, 1 );

  for ( i = 0; i < NumberOfElements; i++ )
    temp(i,0) = vector->data[i * vector->stride];

  return( temp );
}

void Vector_d::SetElements( double Value )
{
  int i;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = Value;
    }
}

double Vector_d::Sum() const
{
  int i;
  double sum = 0.0;

  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];
  return( sum );
}

double Vector_d::Mean() const
{
  int i;
  double sum = 0.0;
  double mean;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];

  mean = sum / NumberOfElements;

  return( mean );
}  

double Vector_d::MaximumElement()
{
  int i;
  double Maximum = -DBL_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > Maximum )
      Maximum = vector->data[i * vector->stride];

  return( Maximum );
}

void Vector_d::RemoveElement( int ElementNumber )
{
  Vector_d temp;

  if ( ElementNumber >= NumberOfElements )
    {
      cout << "WARNING: Vector::RemoveElement Element Number";
      cout << ElementNumber << " > SizeOfArray" << NumberOfElements << " accessed." << endl;
    }

  if ( ElementNumber < 0)
    {
      cout << "WARNING: Vector::RemoveElement ElementNumber";
      cout << ElementNumber << " < 0 accessed." << endl;
    }

  temp.SetNumberOfElements( NumberOfElements );
  for( int i = 0; i < NumberOfElements; i++ )
    temp(i) = gsl_vector_get( vector, i );

  gsl_vector_free( vector );
  NumberOfElements--;
  vector = gsl_vector_calloc( NumberOfElements );
  int counter = -1;
  for ( int i = 0; i < NumberOfElements + 1; i++ ){
    if( i != ElementNumber ){
      counter++;
      gsl_vector_set( vector, counter, temp( i ) );
    }
  }   
}

void Vector_d::AddElement( int ElementNumber )
{
  Vector_d temp;

  if( ElementNumber > NumberOfElements )
    {
      cout << "WARNING: Vector::AddElement Element Number";
      cout << ElementNumber << " > SizeOfArray" << NumberOfElements << " accessed." << endl;
    }

  if( ElementNumber < 0 )
    {
      cout << "WARNING: Vector::AddElement ElementNumber";
      cout << ElementNumber << " < 0 accessed." << endl;
    }

  temp.SetNumberOfElements( NumberOfElements );
  for( int i = 0; i < NumberOfElements; i++ )
    temp(i) = gsl_vector_get( vector, i );

  gsl_vector_free( vector );
  NumberOfElements++;
  vector = gsl_vector_calloc( NumberOfElements );
  int counter = 0;
  for ( int i = 0; i < NumberOfElements - 1; i++ ){
    if( i == ElementNumber ){
      gsl_vector_set( vector, counter, 0 );
      counter++;
    }
    gsl_vector_set( vector, counter, temp( i ) );
    counter++;
  }   
}

