#include "vector.h"

#include "rand.h"
#include "matrix.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "vector_i.h"
#include "matrix_i.h"

#define BOUNDS_CHECK 1

using namespace std;

Vector::Vector():
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
}

Vector::Vector(const Vector& V):
  NumberOfElements(V.NumberOfElements),
  NumberOfMissingValues(V.NumberOfMissingValues),
  MissingValues(0),
  vector(gsl_vector_float_calloc(V.NumberOfElements))
{
  // Copy vector elements
  gsl_vector_float_memcpy(vector, V.vector);

  // Copy missing values
  if (NumberOfMissingValues){
    MissingValues = new int[NumberOfMissingValues];
    int i;
    for (i=0; i<NumberOfMissingValues; i++){
      MissingValues[i] = V.MissingValues[i];
    }
  }
}

Vector::Vector(int NewNumberOfElements):
  NumberOfElements(NewNumberOfElements),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  //? should be an assertion, as
  // new number of elements should
  // probably never be zero
  if ( NewNumberOfElements <= 0 ){
    NumberOfElements = 1;
  }

  vector = gsl_vector_float_calloc( NumberOfElements );
}

Vector::Vector(const char* FileName):
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  Load(FileName);
}

Vector::~Vector()
{
  if (vector != 0){
    gsl_vector_float_free( vector );
  }
  if (NumberOfMissingValues){
    delete [] MissingValues;
  }
}

Vector_d Vector::Double()
{
  Vector_d d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (double)(gsl_vector_float_get( vector, i ) );
   
  return( d );
}
                      
Vector_i Vector::Integer()
{
  Vector_i d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (int)(gsl_vector_float_get( vector, i ) );
   
  return( d );
}
                      
float &Vector::operator()( int element )
{
#ifdef BOUNDS_CHECK
   assert( element < NumberOfElements );
   assert( element > -1 );
#endif
  return( vector->data[element*vector->stride] );
}

float& Vector::operator()(int i) const{
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

float& Vector::operator[](int i){
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

float& Vector::operator[](int i) const{
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

int Vector::GetNumberOfElements() const
{
  return( NumberOfElements );
}

void Vector::SetNumberOfElements( int NewValue )
{
  if ( NumberOfElements == NewValue )
    {
      return;
    }

  //   if ( NumberOfElements != 0 )
  //   {
  if ( vector != 0 )
    gsl_vector_float_free( vector );
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
  vector = gsl_vector_float_calloc( NumberOfElements );
}

FLOAT_ACC Vector::GetElement( int index ) const
{
#ifdef BOUNDS_CHECK   
  if ( index > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector::GetElement.\n";
      exit(0);
    }
  if ( index < 0 )
    {
      cout << "index < 0 in Vector::GetElement.\n";
      exit(0);
    }
#endif
  return( (float)gsl_vector_float_get( vector, index ) );
}

void Vector::SetElement( int index, FLOAT_ACC NewValue )
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
  gsl_vector_float_set( vector, index, NewValue );
}

ostream& operator<< ( ostream& os, const Vector& V )
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

istream& operator>> ( istream& is, Vector& V )
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
	  V( i ) = 0.0;
	  TempMissingValues[V.NumberOfMissingValues++] = i;
	}
      else
	V( i ) = strtod( token, &TailPointer );
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

Vector Vector::operator+( const Vector &v ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to add vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + v.vector->data[offset];
    }

  return( temp );
}

Vector Vector::operator+( double f ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + f;
    }

  return( temp );
}

Vector operator+( double f, const Vector &v )
{
  int i;
  int offset;
  Vector temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = v.vector->data[offset] + f;
    }

  return( temp );
}

Vector &Vector::operator+=( double f )
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

Vector &Vector::operator+=( const Vector &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to add vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] += v.vector->data[offset];
    }

  return( *this );
}

Vector Vector::operator-( const Vector &v ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

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

Vector Vector::operator-( double f ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] - f;
    }

  return( temp );
}

Vector operator-( double f, const Vector &v )
{
  int i;
  int offset;
  Vector temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = f - v.vector->data[offset];
    }

  return( temp );
}

Vector &Vector::operator-=( double f )
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

Vector &Vector::operator-=( const Vector &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to add vectors of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] -= v.vector->data[offset];
    }

  return( *this );
}

Vector Vector::operator-() const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = -vector->data[offset];
    }

  return( temp );
}

Vector Vector::operator*( const Vector &v ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

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

Vector Vector::operator*( double f ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] * f;
    }

  return( temp );
}

Vector operator*( double f, const Vector &v )
{
  int i;
  int offset;
  Vector temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = v.vector->data[offset] * f;
    }

  return( temp );
}

Vector &Vector::operator*=( double f )
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

Vector &Vector::operator*=( const Vector &v )
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

Vector Vector::operator/( const Vector &v ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

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

Vector Vector::operator/( double f ) const
{
  int i;
  int offset;
  Vector temp( NumberOfElements );

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

Vector operator/( double f, const Vector &v )
{
  int i;
  int offset;
  Vector temp( v.NumberOfElements );

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

Vector &Vector::operator/=( double f )
{
  int i;

  if ( f == 0.0 )
    {
      cout << "ERROR: Attempted division of vector by zero." << endl;
      exit ( 1 );
    }
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] /= f;
    }

  return( *this );
}

Vector &Vector::operator/=( const Vector &v )
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

Vector &Vector::operator=( const Vector &v )
{
  int i;

  // Check for self-assignment
   
  if ( &v == this )
    return *this;

  // De-allocate old vector, allocate new one
   
  if ( vector != 0 )
    gsl_vector_float_free( vector );
  NumberOfElements = v.GetNumberOfElements();
  vector = gsl_vector_float_calloc( NumberOfElements );

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
  gsl_vector_float_memcpy( vector, v.vector );

  return *this;
}

int Vector::operator==( const Vector &v ) const
{
  int i;

  if ( NumberOfElements != v.GetNumberOfElements() )
    return( 0 );
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] != v.vector->data[i * v.vector->stride] )
      return( 0 );

  return( 1 );
}

int Vector::GetNumberOfMissingValues() const
{
  return( NumberOfMissingValues );
}

int Vector::GetMissingValueIndex( int MissingValueNumber ) const
{
#ifdef BOUNDS_CHECK   
  if ( MissingValueNumber > NumberOfMissingValues - 1 )
    MissingValueNumber = NumberOfMissingValues - 1;
  if ( MissingValueNumber < 0 )
    MissingValueNumber = 0;
#endif

  return( MissingValues[MissingValueNumber] );
}

int Vector::IsMissingElement( int index ) const
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

void Vector::Load( const char *FileName )
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
    gsl_vector_float_free( vector );
  if ( NumberOfMissingValues != 0 )
    delete [] MissingValues;
      
  NumberOfElements = 0;
  //   }

  // Try to load vector
  
  is >> *this;

  is.close();
}

Matrix Vector::MyTransposeTimesMe() const
{
  int i;
  int j;
  Matrix temp( NumberOfElements, NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      for ( j = 0; j < NumberOfElements; j++ )
	{
	  temp(i,j) = vector->data[i * vector->stride] * vector->data[j * vector->stride];
	}
    }

  return( temp );
}

Matrix Vector::RowMatrix() const
{
  int i;
  Matrix temp( 1, NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    temp(0,i) = vector->data[i * vector->stride];

  return( temp );
}

Matrix Vector::ColumnMatrix() const
{
  int i;
  Matrix temp( NumberOfElements, 1 );

  for ( i = 0; i < NumberOfElements; i++ )
    temp(i,0) = vector->data[i * vector->stride];

  return( temp );
}

void Vector::Randomize( FLOAT_ACC LowerLimit, FLOAT_ACC UpperLimit )
{
  int i;
  double Range = UpperLimit - LowerLimit;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = LowerLimit + Range * myrand();
    }
}

void Vector::RandomizeGaussian( float mu, float sigma )
{
  int i;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = gennor( mu, (double)sigma );
    }
}

void Vector::Sigmoid(void)
{
  int i;
  float overflow=MAX_LIN_ACTIVITY;
  float underflow=MIN_LIN_ACTIVITY;

  for (i = 0; i < NumberOfElements; i++)
    if ( vector->data[i * vector->stride] > overflow )
      vector->data[i * vector->stride] = 1.0;
    else if ( vector->data[i * vector->stride] < underflow )
      vector->data[i * vector->stride] = 0.0;
    else
      vector->data[i * vector->stride] =
	1.0 / ( 1.0 + exp( -vector->data[i * vector->stride] ) );

  return;
}

void Vector::Exp(void)
{
  int i;
  float overflow=MAX_LIN_ACTIVITY;
  float underflow=MIN_LIN_ACTIVITY;

  for (i = 0; i < NumberOfElements; i++)
    if ( vector->data[i * vector->stride] > overflow )
      vector->data[i * vector->stride] = FLT_MAX;
    else if ( vector->data[i * vector->stride] < underflow )
      vector->data[i * vector->stride] = 0.0;
    else
      vector->data[i * vector->stride] =
	exp(vector->data[i * vector->stride]);

  return;
}

void Vector::Tanh(void)
{
  int i;
  float overflow=MAX_LIN_ACTIVITY;
  float underflow=MIN_LIN_ACTIVITY;

  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > overflow )
      vector->data[i * vector->stride] = 1.0;
    else if ( vector->data[i * vector->stride] < underflow )
      vector->data[i * vector->stride] = -1.0;
    else
      vector->data[i * vector->stride] =
	tanh(vector->data[i * vector->stride]);

  return;
}

Vector Vector::Absolute() const
{
  int i;
  Vector temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    temp.SetElement( i, fabs(vector->data[i * vector->stride]));
   
  return( temp );
}

void Vector::SetElements( FLOAT_ACC Value )
{
  int i;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = Value;
    }
}

double Vector::Sum() const
{
  int i;
  double sum = 0.0;

  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];
  return( sum );
}

double Vector::Mean() const
{
  int i;
  double sum = 0.0;
  double mean;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];

  mean = sum / NumberOfElements;

  return( mean );
}  

double Vector::Variance(void) const
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

float Vector::Norm()
{
  int i;
  double sum = 0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride] *
      vector->data[i * vector->stride];
  return( (float)sum );
}

// Functions computing Maximum and Minimum of the elements
float Vector::MinimumElement()
{
  int i;
  float Minimum = FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] < Minimum )
      Minimum = vector->data[i * vector->stride];

  return( Minimum );
}

float Vector::MaximumElement()
{
  int i;
  float Maximum = -FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > Maximum )
      Maximum = vector->data[i * vector->stride];

  return( Maximum );
}

int Vector::WhereMinimumElement()
{
  int i;
  int WhereMinimum = 0;
  float Minimum = FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ ){
    if ( vector->data[i * vector->stride] < Minimum ){
      Minimum = vector->data[i * vector->stride];
      WhereMinimum= i;
    }
  }
  return( WhereMinimum );
}

int Vector::WhereMaximumElement()
{
  int i;
  int WhereMaximum = 0;
  float Maximum = -FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ ){
    if ( vector->data[i * vector->stride] > Maximum ){
      Maximum = vector->data[i * vector->stride];
      WhereMaximum = i;
    }
  }
  return( WhereMaximum );
}

float Vector::MinimumAbsoluteElement()
{
  int i;
  float Minimum = FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( fabs( vector->data[i * vector->stride] ) < Minimum )
      Minimum = fabs( vector->data[i * vector->stride] );

  return( Minimum );
}

float Vector::MaximumAbsoluteElement()
{
  int i;
  float Maximum = 0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( fabs( vector->data[i * vector->stride] ) > Maximum )
      Maximum = fabs( vector->data[i * vector->stride] );

  return( Maximum );
}

int Vector::CumulativeSumDraw() const
{
  int i;
  int SelectedIndex;
  float *Ranges;
  float RandomFloat;
  double TotalSum;
  double CumulativeSum;

  // Calculate sum of elements

  TotalSum = 0.0;
  for ( i = 0; i < NumberOfElements; i++ )
    TotalSum += vector->data[i * vector->stride];

  // Handle pathological case where sum of elements is zero

  if ( TotalSum == 0.0F )
    {
      cout << "Vector::CumulativeSumDraw() ERROR:";
      cout << "Cannot draw from vector where sum of elements is zero." << endl;
      exit( 1 );
    }

  // Ranges[] holds the cumulative sum per element so we can generate
  // a random number in the range 0.0..1.0 and find which bin it
  // falls into. Bins are bigger for larger elements in the array
   
  Ranges = new float[NumberOfElements];

  // Calculate cumulative sum and store scaled version in Ranges[]

  CumulativeSum = 0.0;
  for ( i = 0; i < NumberOfElements; i++ )
    {
      CumulativeSum += vector->data[i * vector->stride] / TotalSum;
      Ranges[i] = CumulativeSum;
    }

  // Random float for bin selection

  RandomFloat = myrand();

  // SelectedIndex is the index corresponding to the Range[] bin that
  // brackets RandomFloat
   
  SelectedIndex = 0;
  for ( i = 1; i < NumberOfElements; i++ )
    if ( RandomFloat > Ranges[i-1] && RandomFloat <= Ranges[i] )
      SelectedIndex = i;

  delete [] Ranges;

  return( SelectedIndex );
}

void Vector::Distinct()
{
  int flag;
  Vector temp( NumberOfElements );

  temp.vector->data[0] = vector->data[0];
  int count = 1;
  for( int i = 1; i < NumberOfElements; i++ ){
    flag = 0;
    for( int j = 0; j < count; j++ )
      if( vector->data[i] == temp.vector->data[j] )
	flag = 1;
    if( flag == 0 ){
      temp.vector->data[count] = vector->data[i];
      count++;}
  }
  NumberOfElements = count;
  for( int i = 0; i < count; i++ )
    vector->data[i] = temp.vector->data[i];
}

void Vector::RemoveElement( int ElementNumber )
{
  Vector temp;

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
    temp(i) = gsl_vector_float_get( vector, i );

  gsl_vector_float_free( vector );
  NumberOfElements--;
  vector = gsl_vector_float_calloc( NumberOfElements );
  int counter = -1;
  for ( int i = 0; i < NumberOfElements + 1; i++ ){
    if( i != ElementNumber ){
      counter++;
      gsl_vector_float_set( vector, counter, temp( i ) );
    }
  }   
}

void Vector::AddElement( int ElementNumber )
{
  Vector temp;

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
    temp(i) = gsl_vector_float_get( vector, i );

  gsl_vector_float_free( vector );
  NumberOfElements++;
  vector = gsl_vector_float_calloc( NumberOfElements );
  int counter = 0;
  for ( int i = 0; i < NumberOfElements - 1; i++ ){
    if( i == ElementNumber ){
      gsl_vector_float_set( vector, counter, 0 );
      counter++;
    }
    gsl_vector_float_set( vector, counter, temp( i ) );
    counter++;
  }   
}

void Vector::Sort()
{
  float v[ NumberOfElements ];

  for( int j = 0; j < NumberOfElements; j++ )
    v[j] = gsl_vector_float_get( vector, j );

  qsort( v, NumberOfElements, sizeof (float), CompareVector );

  for( int j = 0; j < NumberOfElements; j++ )
    gsl_vector_float_set( vector, j, v[j] );
}

int CompareVector( const void *cv1, const void *cv2 )
{
  const float *key1 = (const float *)cv1;
  const float *key2 = (const float *)cv2;
   
  if ( *key1 < *key2 )
    return( -1 );
  else if ( *key2 < *key1 )
    return( 1 );
  else
    return( 0 );
}

gsl_vector_float* Vector::GetGslVector()
{
  return( vector );
}
