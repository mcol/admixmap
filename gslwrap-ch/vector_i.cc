#include "vector_i.h"

#include "rand.h"
#include "vector.h"
#include "matrix.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "matrix_i.h"

#define BOUNDS_CHECK 1

using namespace std;;

Vector_i::Vector_i():
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
}

Vector_i::Vector_i(const Vector_i& V):
  NumberOfElements(V.NumberOfElements),
  NumberOfMissingValues(V.NumberOfMissingValues),
  MissingValues(0),
  vector(gsl_vector_int_calloc(V.NumberOfElements))
{

  // Copy vector elements
  gsl_vector_int_memcpy(vector,V.vector);

  // Copy missing values
  if (NumberOfMissingValues){
    int i;
    MissingValues = new int[NumberOfMissingValues];
    for (i=0; i<NumberOfMissingValues; i++){
      MissingValues[i] = V.GetMissingValueIndex( i );
    }
  }
}

Vector_i::Vector_i(int NewNumberOfElements):
  NumberOfElements(NewNumberOfElements),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  //? should be an assertion
  if (NewNumberOfElements<=0){
    NumberOfElements = 1;
  }

  vector = gsl_vector_int_calloc(NumberOfElements);
}

Vector_i::Vector_i(const char* FileName):
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  Load( FileName );
}

Vector_i::~Vector_i()
{
  if ( vector != 0 )
    gsl_vector_int_free( vector );
  if ( NumberOfMissingValues != 0 )
    delete [] MissingValues;
}

Vector Vector_i::Float()
{
  Vector d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (float)(gsl_vector_int_get( vector, i ) );
   
  return( d );
}

Vector_d Vector_i::Double()
{
  Vector_d d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (double)(gsl_vector_int_get( vector, i ) );
   
  return( d );
}
                      
int &Vector_i::operator()( int element )
{
#ifdef BOUNDS_CHECK   
  if ( element > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector_i::().\n";
      exit(0);
    }
  if ( element < 0 )
    {
      cout << "index < 0 in Vector_i::().\n";
      exit(0);
    }
#endif
  return( vector->data[element*vector->stride] );
}

int& Vector_i::operator()(int i) const{
#ifdef BOUNDS_CHECK   
  if ( i > NumberOfElements - 1 ){
    cout << "index = " << i << " > NumberOfElements - 1 in Vector_i::().\n";
    exit(0);
  }
  if ( i < 0 ){
    cout << "index = " << i << " < 0 in Vector_i::().\n";
    exit(0);
  }
#endif
  return( vector->data[i*vector->stride] );  
}

int& Vector_i::operator[](int i){
#ifdef BOUNDS_CHECK   
  if ( i > NumberOfElements - 1 ){
    cout << "index = " << i << " > NumberOfElements - 1 in Vector_i::[].\n";
    exit(0);
  }
  if ( i < 0 ){
    cout << "index = " << i << " < 0 in Vector_i::[].\n";
    exit(0);
  }
#endif
  return( vector->data[i*vector->stride] );  
}

int& Vector_i::operator[](int i) const{
#ifdef BOUNDS_CHECK   
  if ( i > NumberOfElements - 1 ){
    cout << "index = " << i << " > NumberOfElements - 1 in Vector_i::[].\n";
    exit(0);
  }
  if ( i < 0 ){
    cout << "index = " << i << " < 0 in Vector_i::[].\n";
    exit(0);
  }
#endif
  return( vector->data[i*vector->stride] );  
}

int Vector_i::GetNumberOfElements() const
{
  return( NumberOfElements );
}

void Vector_i::SetNumberOfElements( int NewValue )
{
  if ( NumberOfElements == NewValue )
    {
      return;
    }

  //   if ( NumberOfElements != 0 )
  //   {
  if ( vector != 0 )
    gsl_vector_int_free( vector );
  //   }

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
  vector = gsl_vector_int_calloc( NumberOfElements );
}

int Vector_i::GetElement( int index ) const
{
#ifdef BOUNDS_CHECK   
  if ( index > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector_i::().\n";
      exit(0);
    }
  if ( index < 0 )
    {
      cout << "index < 0 in Vector_i::().\n";
      exit(0);
    }
#endif
  return( (gsl_vector_int_get( vector, index ) ) );
}

void Vector_i::SetElement( int index, int NewValue )
{
#ifdef BOUNDS_CHECK   
  if ( NumberOfElements == 0 )
    {
      return;
    }
  if ( index > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector::SetElement.\n";
      exit(0);
    }
  if ( index < 0 )
    {
      cout << "index < 0 in Vector::SetElement.\n";
      exit(0);
    }
#endif
  gsl_vector_int_set( vector, index, NewValue );
}

ostream& operator<< ( ostream& os, const Vector_i& V )
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

istream& operator>> ( istream& is, Vector_i& V )
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
      //       cout << "CurrentLine:" << endl << CurrentLine << endl;
      //       cout << "Token: \"" << token << "\"" << endl;
      if ( !strcmp( token, "#" ) )
	{
	  V( i ) = 0;
	  TempMissingValues[V.NumberOfMissingValues++] = i;
	}
      else
	V( i ) = (int)strtod( token, &TailPointer );
      //       cout << "Double conversion of token: " << V( i ) << endl;
      //       cout << "Missing values: " << V.NumberOfMissingValues << endl;
      i++;
    }

  // Copy missing value indices from TempMissingValues to V.MissingValues
   
  if ( V.NumberOfMissingValues )
    V.MissingValues = new int[V.NumberOfMissingValues];
  else
    V.MissingValues = 0;
  for ( i = 0; i < V.NumberOfMissingValues; i++ )
    V.MissingValues[i] = TempMissingValues[i];
  //    for ( i = 0; i < V.NumberOfMissingValues; i++ )
  //       cout << "Missing value[" << i << "] = " << V.MissingValues[i] << endl;

  // Peek one character ahead to test for eof, put char back if we're not at eof
   
  is.get( PreviousChar );
  if ( is )
    is.putback( PreviousChar );
   
  delete [] CurrentLine;
  delete [] TempMissingValues;
   
  return is;
}	

Vector_i Vector_i::operator+( const Vector_i &v ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR in vector_i: Attempt to add vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + v.vector->data[offset];
    }

  return( temp );
}

Vector_i Vector_i::operator+( int f ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + f;
    }

  return( temp );
}

Vector_i operator+( int f, const Vector_i &v )
{
  int i;
  int offset;
  Vector_i temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = v.vector->data[offset] + f;
    }

  return( temp );
}

Vector_i &Vector_i::operator+=( int f )
{
  int i;
  int offset;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] += f;
    }

  return( *this );
}

Vector_i &Vector_i::operator+=( const Vector_i &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR in vector_i: Attempt to add vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] += v.vector->data[offset];
    }

  return( *this );
}

Vector_i Vector_i::operator-( const Vector_i &v ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to subtract vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] - v.vector->data[offset];
    }

  return( temp );
}

Vector_i Vector_i::operator-( int f ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] - f;
    }

  return( temp );
}

Vector_i operator-( int f, const Vector_i &v )
{
  int i;
  int offset;
  Vector_i temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = f - v.vector->data[offset];
    }

  return( temp );
}

Vector_i &Vector_i::operator-=( int f )
{
  int i;
  int offset;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] -= f;
    }

  return( *this );
}

Vector_i &Vector_i::operator-=( const Vector_i &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR in vector_i: Attempt to add vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] -= v.vector->data[offset];
    }

  return( *this );
}

Vector_i Vector_i::operator-() const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = -vector->data[offset];
    }

  return( temp );
}

Vector_i Vector_i::operator*( const Vector_i &v ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to multiply vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] * v.vector->data[offset];
    }

  return( temp );
}

Vector_i Vector_i::operator*( int f ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] * f;
    }

  return( temp );
}

Vector_i operator*( int f, const Vector_i &v )
{
  int i;
  int offset;
  Vector_i temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = v.vector->data[offset] * f;
    }

  return( temp );
}

Vector_i &Vector_i::operator*=( int f )
{
  int i;
  int offset;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] *= f;
    }

  return( *this );
}

Vector_i &Vector_i::operator*=( const Vector_i &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to multiply vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] *= v.vector->data[offset];
    }

  return( *this );
}

Vector_i Vector_i::operator/( const Vector_i &v ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to divide vectors of different length." << endl;
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

Vector_i Vector_i::operator/( int f ) const
{
  int i;
  int offset;
  Vector_i temp( NumberOfElements );

  if ( f == 0.0 )
    {
      cout << "ERROR: Attempted to divide vector by zero." << endl;
      exit( 1 );
    }
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] / f;
    }

  return( temp );
}

Vector_i operator/( int f, const Vector_i &v )
{
  int i;
  int offset;
  Vector_i temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      if ( v.vector->data[offset] == 0.0 )
	{
	  cout << "ERROR: Attempted to divide by zero element of a vector.";
	  cout << endl;
	  exit( 1 );
	}
      temp.vector->data[offset] = f / v.vector->data[offset];
    }

  return( temp );
}

Vector_i &Vector_i::operator/=( int f )
{
  int i;

  if ( f == 0.0 )
    {
      cout << "ERROR: Attempted division of vector_i by zero." << endl;
      exit ( 1 );
    }
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] /= f;
    }

  return( *this );
}

Vector_i &Vector_i::operator/=( const Vector_i &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to divide vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      if ( v.vector->data[offset] == 0.0 )
	{
	  cout << "ERROR: Attempted to divide by zero element of a vector." << endl;
	  exit( 1 );
	}
      vector->data[offset] /= v.vector->data[offset];
    }

  return( *this );
}

Vector_i &Vector_i::operator=( const Vector_i &v )
{
  int i;

  // Check for self-assignment
   
  if ( &v == this )
    return *this;

  // De-allocate old vector, allocate new one
   
  if ( vector != 0 )
    gsl_vector_int_free( vector );
  NumberOfElements = v.GetNumberOfElements();
  vector = gsl_vector_int_calloc( NumberOfElements );

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
  gsl_vector_int_memcpy( vector, v.vector );

  return *this;
}

int Vector_i::operator==( const Vector_i &v ) const
{
  int i;

  if ( NumberOfElements != v.GetNumberOfElements() )
    return( 0 );
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] != v.vector->data[i * v.vector->stride] )
      return( 0 );

  return( 1 );
}

int Vector_i::GetNumberOfMissingValues() const
{
  return( NumberOfMissingValues );
}

int Vector_i::GetMissingValueIndex( int MissingValueNumber ) const
{
#ifdef BOUNDS_CHECK   
  if ( MissingValueNumber > NumberOfMissingValues - 1 )
    MissingValueNumber = NumberOfMissingValues - 1;
  if ( MissingValueNumber < 0 )
    MissingValueNumber = 0;
#endif

  return( MissingValues[MissingValueNumber] );
}

int Vector_i::IsMissingElement( int index ) const
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

void Vector_i::Load( const char *FileName )
{
  ifstream is( FileName );

  // Check that file exists
  
  if ( !is )
    {
      cout << "ERROR: Couldn't open file " << FileName << endl;
    }

  // If vector already allocated, deallocate

  //   if ( NumberOfElements != 0 )
  //   {
  if ( vector != 0 )
    gsl_vector_int_free( vector );
  if ( NumberOfMissingValues != 0 )
    delete [] MissingValues;
      
  NumberOfElements = 0;
  //   }

  // Try to load vector
  
  is >> *this;

  is.close();
}

Matrix_i Vector_i::MyTransposeTimesMe() const
{
  int i;
  int j;
  Matrix_i temp( NumberOfElements, NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      for ( j = 0; j < NumberOfElements; j++ )
	{
	  temp(i,j) = vector->data[i * vector->stride] * vector->data[j * vector->stride];
	}
    }

  return( temp );
}

Matrix_i Vector_i::RowMatrix() const
{
  int i;
  Matrix_i temp( 1, NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    temp(0,i) = vector->data[i * vector->stride];

  return( temp );
}

Matrix_i Vector_i::ColumnMatrix() const
{
  int i;
  Matrix_i temp( NumberOfElements, 1 );

  for ( i = 0; i < NumberOfElements; i++ )
    temp(i,0) = vector->data[i * vector->stride];

  return( temp );
}

void Vector_i::Randomize( int LowerLimit, int UpperLimit )
{
  int i;
  int Range = UpperLimit - LowerLimit + 1;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = LowerLimit + Range * (int)myrand();
    }
}


Vector_i Vector_i::Absolute() const
{
  int i;
  Vector_i temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ ){
    temp.SetElement( i, static_cast<int>(fabs((double)vector->data[i * vector->stride])));
  }
  return( temp );
}

void Vector_i::SetElements( int Value )
{
  int i;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = Value;
    }
}

int Vector_i::Sum() const
{
  int i;
  int sum = 0;

  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];
  return( sum );
}

double Vector_i::Mean() const
{
  int i;
  double sum = 0.0;
  double mean;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];

  mean = sum / NumberOfElements;

  return( mean );
}  

double Vector_i::Variance(void) const
{
  int i;
  double mean = Mean();
  double variance = 0.0;
  double app = 0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      double s = vector->data[i * vector->stride] - mean;
      variance += s * s;
      app += s;
    }

  return((variance - (app * app) / NumberOfElements) / (NumberOfElements - 1));
}

double Vector_i::Norm()
{
  int i;
  double sum = 0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride] *
      vector->data[i * vector->stride];
  return( (double)sum );
}

// Functions computing Maximum and Minimum of the elements
int Vector_i::MinimumElement()
{
  int i;
  int Minimum = INT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] < Minimum )
      Minimum = vector->data[i * vector->stride];

  return( Minimum );
}

int Vector_i::MaximumElement()
{
  int i;
  int Maximum = -INT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > Maximum )
      Maximum = vector->data[i * vector->stride];

  return( Maximum );
}

int Vector_i::MinimumAbsoluteElement()
{
  int i;
  double Minimum = INT_MAX;
  double result = 0.0;

  for ( i = 0; i < NumberOfElements; i++ ){
    result = fabs((double)vector->data[i * vector->stride] );
    if ( result < Minimum ){
      Minimum = result;
    }
  }
  return( static_cast<int>(Minimum) );
}

int Vector_i::MaximumAbsoluteElement()
{
  int i;
  double Maximum = 0.0;
  double result = 0.0;
  
  for ( i = 0; i < NumberOfElements; i++ ){
    result = fabs((double)vector->data[i * vector->stride] );
    if ( result > Maximum ){
      Maximum = result;
    }
  }
  return( static_cast<int>(Maximum) );
}

void Vector_i::Distinct()
{
  int flag;
  Vector_i temp( NumberOfElements );

  temp.vector->data[0] = vector->data[0];
  int count = 1;
  for( int i = 1; i < NumberOfElements; i++ ){
    flag = 0;
    for( int j = 0; j < count; j++ )
      if( vector->data[i] == temp.vector->data[j] )
	flag = 1;
    if( flag == 0 ){
      temp.vector->data[count] = vector->data[i];
      count++;
    }
  }
  gsl_vector_int_free( vector );
  NumberOfElements = count;
  vector = gsl_vector_int_calloc( NumberOfElements );
  for( int i = 0; i < count; i++ )
    vector->data[i] = temp.vector->data[i];
}

void Vector_i::AddElement( int ElementNumber )
{
  Vector_i temp;

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
    temp(i) = gsl_vector_int_get( vector, i );

  gsl_vector_int_free( vector );
  NumberOfElements++;
  vector = gsl_vector_int_calloc( NumberOfElements );
  int counter = 0;
  for ( int i = 0; i < NumberOfElements - 1; i++ ){
    if( i == ElementNumber ){
      gsl_vector_int_set( vector, counter, 0 );
      counter++;
    }
    gsl_vector_int_set( vector, counter, temp( i ) );
    counter++;
  }   
}

void Vector_i::RemoveElement( int ElementNumber )
{
  Vector_i temp;

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
    temp(i) = gsl_vector_int_get( vector, i );

  gsl_vector_int_free( vector );
  NumberOfElements--;
  vector = gsl_vector_int_calloc( NumberOfElements );
  int counter = -1;
  for ( int i = 0; i < NumberOfElements + 1; i++ ){
    if( i != ElementNumber ){
      counter++;
      gsl_vector_int_set( vector, counter, temp( i ) );
    }
  }   
}

gsl_vector_int* Vector_i::GetGslVector()
{
  return( vector );
}
