#include "rand.h"

#include "vector.h"
#include "vector_d.h"
#include "vector_i.h"

#include "matrix.h"
#include "matrix_d.h"
#include "matrix_i.h"

#define BOUNDS_CHECK 1

using namespace std;

Vector%TYPE2%::Vector%TYPE2%():
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
}

Vector%TYPE2%::Vector%TYPE2%(const Vector%TYPE2%& V):
  NumberOfElements(V.NumberOfElements),
  NumberOfMissingValues(V.NumberOfMissingValues),
  MissingValues(0),
  vector(gsl_vector%GSL%_calloc(V.NumberOfElements))
{

  // Copy vector elements
  gsl_vector%GSL%_memcpy( vector, V.vector );

  // Copy missing values
  if (NumberOfMissingValues){
    int i;
    MissingValues = new int[NumberOfMissingValues];
    for (i=0; i<NumberOfMissingValues; i++){
      MissingValues[i] = V.MissingValues[i];
    }
  }
}

Vector%TYPE2%::Vector%TYPE2%(int NewNumberOfElements):
  NumberOfElements(NewNumberOfElements),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  //? should be an assertion
  if ( NewNumberOfElements <= 0 ){
    NumberOfElements = 1;
  }

  vector = gsl_vector%GSL%_calloc( NumberOfElements );
}

Vector%TYPE2%::Vector%TYPE2%(const char* FileName):
  NumberOfElements(0),
  NumberOfMissingValues(0),
  MissingValues(0),
  vector(0)
{
  Load(FileName);
}

Vector%TYPE2%::~Vector%TYPE2%()
{
  if (vector != 0){
    gsl_vector%GSL%_free(vector);
  }
  if (NumberOfMissingValues){
    delete [] MissingValues;
  }
}

Vector Vector%TYPE2%::Float()
{
  Vector d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (float)(gsl_vector%GSL%_get( vector, i ) );
   
  return( d );
}
                      
Vector_d Vector%TYPE2%::Double()
{
  Vector_d d( NumberOfElements );

  for( int i = 0; i < NumberOfElements; i++ )
    d(i) = (double)(gsl_vector%GSL%_get( vector, i ) );
   
  return( d );
}
                      
%TYPE% &Vector%TYPE2%::operator()( int element )
{
#ifdef BOUNDS_CHECK   
  if ( element > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector%TYPE2%::().\n";
      exit(0);
    }
  if ( element < 0 )
    {
      cout << "index < 0 in Vector%TYPE2%::().\n";
      exit(0);
    }
#endif
  return( vector->data[element*vector->stride] );
}

%TYPE%& Vector%TYPE2%::operator()(int i) const{
#ifdef BOUNDS_CHECK   
  if ( i > NumberOfElements - 1 ){
    cout << "index = " << i << " > NumberOfElements - 1 in Vector%TYPE2%::().\n";
    exit(0);
  }
  if ( i < 0 ){
    cout << "index = " << i << " < 0 in Vector_i::().\n";
    exit(0);
  }
#endif
  return( vector->data[i*vector->stride] );  
}

%TYPE%& Vector%TYPE2%::operator[](int i){
#ifdef BOUNDS_CHECK   
  if ( i > NumberOfElements - 1 ){
    cout << "index = " << i << " > NumberOfElements - 1 in Vector%TYPE2%::[].\n";
    exit(0);
  }
  if ( i < 0 ){
    cout << "index = " << i << " < 0 in Vector_i::[].\n";
    exit(0);
  }
#endif
  return( vector->data[i*vector->stride] );  
}

%TYPE%& Vector%TYPE2%::operator[](int i) const{
#ifdef BOUNDS_CHECK   
  if ( i > NumberOfElements - 1 ){
    cout << "index = " << i << " > NumberOfElements - 1 in Vector%TYPE2%::[].\n";
    exit(0);
  }
  if ( i < 0 ){
    cout << "index = " << i << " < 0 in Vector_i::[].\n";
    exit(0);
  }
#endif
  return( vector->data[i*vector->stride] );  
}

int Vector%TYPE2%::GetNumberOfElements() const
{
  return( NumberOfElements );
}

void Vector%TYPE2%::SetNumberOfElements( int NewValue )
{
  if ( NumberOfElements == NewValue )
    {
      return;
    }

  if ( vector != 0 ){
    gsl_vector%GSL%_free( vector );
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
  vector = gsl_vector%GSL%_calloc( NumberOfElements );
}

%TYPE% Vector%TYPE2%::GetElement( int index ) const
{
#ifdef BOUNDS_CHECK   
  if ( index > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector%TYPE2%::GetElement.\n";
      exit(0);
    }
  if ( index < 0 )
    {
      cout << "index < 0 in Vector%TYPE2%::GetElement.\n";
      exit(0);
    }
#endif
  return( (%TYPE%)gsl_vector%GSL%_get( vector, index ) );
}

void Vector%TYPE2%::SetElement( int index, %TYPE% NewValue )
{
#ifdef BOUNDS_CHECK   
  if ( index > NumberOfElements - 1 )
    {
      cout << "index > NumberOfElements - 1 in Vector%TYPE2%::SetElement.\n";
      exit(0);
    }
  if ( index < 0 )
    {
      cout << "index < 0 in Vector%TYPE2%::SetElement.\n";
      exit(0);
    }
  if ( NumberOfElements == 0 )
    {
      return;
    }
#endif
  gsl_vector%GSL%_set( vector, index, NewValue );
}

ostream& operator<< ( ostream& os, const Vector%TYPE2%& V )
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

istream& operator>> ( istream& is, Vector%TYPE2%& V )
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
	V( i ) = (%TYPE%)0.0;
	TempMissingValues[V.NumberOfMissingValues++] = i;
      } else {
	V( i ) = (%TYPE%)strtod( token, &TailPointer );
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

Vector%TYPE2% Vector%TYPE2%::operator+( const Vector%TYPE2% &v ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to add vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + v.vector->data[offset];
    }

  return( temp );
}

Vector%TYPE2% Vector%TYPE2%::operator+( %TYPE% f ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] + f;
    }

  return( temp );
}

Vector%TYPE2% operator+( %TYPE% f, const Vector%TYPE2% &v )
{
  int i;
  int offset;
  Vector%TYPE2% temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = v.vector->data[offset] + f;
    }

  return( temp );
}

Vector%TYPE2% &Vector%TYPE2%::operator+=( %TYPE% f )
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

Vector%TYPE2% &Vector%TYPE2%::operator+=( const Vector%TYPE2% &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to add vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] += v.vector->data[offset];
    }

  return( *this );
}

Vector%TYPE2% Vector%TYPE2%::operator-( const Vector%TYPE2% &v ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

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

Vector%TYPE2% Vector%TYPE2%::operator-( %TYPE% f ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] - f;
    }

  return( temp );
}

Vector%TYPE2% operator-( %TYPE% f, const Vector%TYPE2% &v )
{
  int i;
  int offset;
  Vector%TYPE2% temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = f - v.vector->data[offset];
    }

  return( temp );
}

Vector%TYPE2% &Vector%TYPE2%::operator-=( %TYPE% f )
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

Vector%TYPE2% &Vector%TYPE2%::operator-=( const Vector%TYPE2% &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to add vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] -= v.vector->data[offset];
    }

  return( *this );
}

Vector%TYPE2% Vector%TYPE2%::operator-() const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = -vector->data[offset];
    }

  return( temp );
}

Vector%TYPE2% Vector%TYPE2%::operator*( const Vector%TYPE2% &v ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to multiply vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] * v.vector->data[offset];
    }

  return( temp );
}

Vector%TYPE2% Vector%TYPE2%::operator*( %TYPE% f ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] * f;
    }

  return( temp );
}

Vector%TYPE2% operator*( %TYPE% f, const Vector%TYPE2% &v )
{
  int i;
  int offset;
  Vector%TYPE2% temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      temp.vector->data[offset] = v.vector->data[offset] * f;
    }

  return( temp );
}

Vector%TYPE2% &Vector%TYPE2%::operator*=( %TYPE% f )
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

Vector%TYPE2% &Vector%TYPE2%::operator*=( const Vector%TYPE2% &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to multiply vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      vector->data[offset] *= v.vector->data[offset];
    }

  return( *this );
}

Vector%TYPE2% Vector%TYPE2%::operator/( const Vector%TYPE2% &v ) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to divide vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      if ( v.vector->data[offset] == (%TYPE%)0.0 )
	{
	  cout << "ERROR: Division by zero in element " << i << " of denominator vector." << endl;
	  exit( 1 );
	}
      temp.vector->data[offset] = vector->data[offset] / v.vector->data[offset];
    }

  return( temp );
}

Vector%TYPE2% Vector%TYPE2%::operator/(double f) const
{
  int i;
  int offset;
  Vector%TYPE2% temp( NumberOfElements );

  if ( f == 0.0 )
    {
      cout << "ERROR: Attempted to divide vector_d by zero." << endl;
      exit( 1 );
    }
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * temp.vector->stride;
      temp.vector->data[offset] = vector->data[offset] / f;
    }

  return( temp );
}

Vector%TYPE2% operator/( %TYPE% f, const Vector%TYPE2% &v )
{
  int i;
  int offset;
  Vector%TYPE2% temp( v.NumberOfElements );

  for ( i = 0; i < v.NumberOfElements; i++ )
    {
      offset = i * v.vector->stride;
      if ( v.vector->data[offset] == (%TYPE%)0.0 )
	{
	  cout << "ERROR: Attempted to divide by zero element of a vector.";
	  cout << endl;
	  exit( 1 );
	}
      temp.vector->data[offset] = f / v.vector->data[offset];
    }

  return( temp );
}

Vector%TYPE2% &Vector%TYPE2%::operator/=( %TYPE% f )
{
  int i;

  if ( f == (%TYPE%)0.0 )
    {
      cout << "ERROR: Attempted division of vector_d by zero." << endl;
      exit ( 1 );
    }
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] /= f;
    }

  return( *this );
}

Vector%TYPE2% &Vector%TYPE2%::operator/=( const Vector%TYPE2% &v )
{
  int i;
  int offset;

  if ( v.NumberOfElements != NumberOfElements )
    {
      cout << "ERROR: Attempt to divide vectors(_d) of different length." << endl;
      exit( 1 );
    }

  for ( i = 0; i < NumberOfElements; i++ )
    {
      offset = i * vector->stride;
      if ( v.vector->data[offset] == (%TYPE%)0.0 )
	{
	  cout << "ERROR: Attempted to divide by zero element of a vector(_d)." << endl;
	  exit( 1 );
	}
      vector->data[offset] /= v.vector->data[offset];
    }

  return( *this );
}

Vector%TYPE2% &Vector%TYPE2%::operator=( const Vector%TYPE2% &v )
{
  int i;

  // Check for self-assignment
   
  if ( &v == this )
    return *this;

  // De-allocate old vector, allocate new one
   
  if ( vector != 0 )
    gsl_vector%GSL%_free( vector );
  NumberOfElements = v.GetNumberOfElements();
  vector = gsl_vector%GSL%_calloc( NumberOfElements );

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
  gsl_vector%GSL%_memcpy( vector, v.vector );

  return *this;
}

int Vector%TYPE2%::operator==( const Vector%TYPE2% &v ) const
{
  int i;

  if ( NumberOfElements != v.GetNumberOfElements() )
    return( 0 );
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] != v.vector->data[i * v.vector->stride] )
      return( 0 );

  return( 1 );
}

int Vector%TYPE2%::GetNumberOfMissingValues() const
{
  return( NumberOfMissingValues );
}

int Vector%TYPE2%::GetMissingValueIndex( int MissingValueNumber ) const
{
#ifdef BOUNDS_CHECK   
  if ( MissingValueNumber > NumberOfMissingValues - 1 )
    MissingValueNumber = NumberOfMissingValues - 1;
  if ( MissingValueNumber < 0 )
    MissingValueNumber = 0;
#endif

  return( MissingValues[MissingValueNumber] );
}

int Vector%TYPE2%::IsMissingElement( int index ) const
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

void Vector%TYPE2%::Load( const char *FileName )
{
  ifstream is( FileName );

  // Check that file exists
  
  if ( !is )
    {
      cout << "ERROR: Couldn't open file " << FileName << endl;
    }

  // If vector already allocated, deallocate

  if ( vector != 0 )
    gsl_vector%GSL%_free( vector );
  if ( NumberOfMissingValues != 0 )
    delete [] MissingValues;
      
  NumberOfElements = 0;

  // Try to load vector
  
  is >> *this;

  is.close();
}
/*Next three functions cause compiler errors
Matrix%TYPE2% Vector%TYPE2%::MyTransposeTimesMe() const
{
  int i;
  int j;
  Matrix%TYPE2% temp( NumberOfElements, NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    {
      for ( j = 0; j < NumberOfElements; j++ )
	{
	  temp.matrix->data[i * temp.matrix->size2 + j] = vector->data[i * vector->stride] * vector->data[j * vector->stride];
	}
    }

  return( temp );
}

Matrix%TYPE2% Vector%TYPE2%::RowMatrix() const
{
  int i;
  Matrix%TYPE2% temp( 1, NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    temp.matrix->data[i] = vector->data[i * vector->stride];

  return( temp );
}

Matrix%TYPE2% Vector%TYPE2%::ColumnMatrix() const
{
  int i;
  Matrix%TYPE2% temp( NumberOfElements, 1 );

  for ( i = 0; i < NumberOfElements; i++ )
    temp.matrix->data[i * temp.matrix->size2] = vector->data[i * vector->stride];

  return( temp );
}
*/
void Vector%TYPE2%::Randomize( %TYPE% LowerLimit, %TYPE% UpperLimit )
{
  int i;
  %TYPE% Range = UpperLimit - LowerLimit;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = LowerLimit + Range * myrand();
    }
}

void Vector%TYPE2%::RandomizeGaussian( %TYPE% mu, %TYPE% sigma )
{
  int i;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = (%TYPE%)gennor( mu, (%TYPE%)sigma );
    }
}

void Vector%TYPE2%::Sigmoid(void)
{
  int i;
  %TYPE% overflow  = (%TYPE%)MAX_LIN_ACTIVITY;
  %TYPE% underflow = (%TYPE%)MIN_LIN_ACTIVITY;

  for (i = 0; i < NumberOfElements; i++)
    if ( vector->data[i * vector->stride] > overflow )
      vector->data[i * vector->stride] = (%TYPE%)1.0;
    else if ( vector->data[i * vector->stride] < underflow )
      vector->data[i * vector->stride] = (%TYPE%)0.0;
    else
      vector->data[i * vector->stride] = (%TYPE%)
	1.0 / ( 1.0 + exp((double) -vector->data[i * vector->stride] ) );

  return;
}

void Vector%TYPE2%::Exp(void)
{
  int i;
  %TYPE% overflow  = (%TYPE%)MAX_LIN_ACTIVITY;
  %TYPE% underflow = (%TYPE%)MIN_LIN_ACTIVITY;

  for (i = 0; i < NumberOfElements; i++)
    if ( vector->data[i * vector->stride] > overflow )
      vector->data[i * vector->stride] = FLT_MAX;
    else if ( vector->data[i * vector->stride] < underflow )
      vector->data[i * vector->stride] = (%TYPE%)0.0;
    else
      vector->data[i * vector->stride] =
	exp((double)vector->data[i * vector->stride]);

  return;
}

void Vector%TYPE2%::Tanh(void)
{
  int i;
  %TYPE% overflow  = (%TYPE%)MAX_LIN_ACTIVITY;
  %TYPE% underflow = (%TYPE%)MIN_LIN_ACTIVITY;

  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > overflow )
      vector->data[i * vector->stride] = 1.0;
    else if ( vector->data[i * vector->stride] < underflow )
      vector->data[i * vector->stride] = -1.0;
    else
      vector->data[i * vector->stride] =
	tanh((double)vector->data[i * vector->stride]);

  return;
}

Vector%TYPE2% Vector%TYPE2%::Absolute() const
{
  int i;
  Vector%TYPE2% temp( NumberOfElements );

  for ( i = 0; i < NumberOfElements; i++ )
    temp.SetElement( i, fabs((double)vector->data[i * vector->stride]));
   
  return( temp );
}

void Vector%TYPE2%::SetElements( %TYPE% Value )
{
  int i;

  for ( i = 0; i < NumberOfElements; i++ )
    {
      vector->data[i * vector->stride] = Value;
    }
}

%TYPE% Vector%TYPE2%::Sum() const
{
  int i;
  %TYPE% sum = (%TYPE%)0.0;

  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];
  return( sum );
}

%TYPE% Vector%TYPE2%::Mean() const
{
  int i;
  %TYPE% sum = (%TYPE%)0.0;
  %TYPE% mean;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride];

  mean = sum / NumberOfElements;

  return( mean );
}  

%TYPE% Vector%TYPE2%::Variance(void) const
{
  int i;
  %TYPE% mean = Mean();
  %TYPE% variance = (%TYPE%)0.0;
  %TYPE% app = (%TYPE%)0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    {
      %TYPE% s = vector->data[i * vector->stride] - mean;
      variance += s * s;
      app += s;
    }

  return((variance - (app * app) / NumberOfElements) / (NumberOfElements - 1));
}

%TYPE% Vector%TYPE2%::Norm()
{
  int i;
  %TYPE% sum = (%TYPE%)0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    sum += vector->data[i * vector->stride] *
      vector->data[i * vector->stride];
  return( (%TYPE%)sum );
}

// Functions computing Maximum and Minimum of the elements
%TYPE% Vector%TYPE2%::MinimumElement()
{
  int i;
  %TYPE% Minimum = DBL_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] < Minimum )
      Minimum = vector->data[i * vector->stride];

  return( Minimum );
}

%TYPE% Vector%TYPE2%::MaximumElement()
{
  int i;
  %TYPE% Maximum = -DBL_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > Maximum )
      Maximum = vector->data[i * vector->stride];

  return( Maximum );
}

int Vector%TYPE2%::WhereMinimumElement()
{
  int i;
  int WhereMinimum = 0;
  float Minimum = FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] < Minimum ){
      Minimum = vector->data[i * vector->stride];
      WhereMinimum= i;
    }

  return( WhereMinimum );
}

int Vector%TYPE2%::WhereMaximumElement()
{
  int i;
  int WhereMaximum = 0;
  float Maximum = -FLT_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( vector->data[i * vector->stride] > Maximum ){
      Maximum = vector->data[i * vector->stride];
      WhereMaximum = i;
    }

  return( WhereMaximum );
}

%TYPE% Vector%TYPE2%::MinimumAbsoluteElement()
{
  int i;
  %TYPE% Minimum = DBL_MAX;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( fabs((double)vector->data[i * vector->stride] ) < Minimum )
      Minimum = fabs((double)vector->data[i * vector->stride] );

  return( Minimum );
}

%TYPE% Vector%TYPE2%::MaximumAbsoluteElement()
{
  int i;
  %TYPE% Maximum = (%TYPE%)0.0;
   
  for ( i = 0; i < NumberOfElements; i++ )
    if ( fabs((double)vector->data[i * vector->stride] ) > Maximum )
      Maximum = fabs((double)vector->data[i * vector->stride] );

  return( Maximum );
}

int Vector%TYPE2%::CumulativeSumDraw() const
{
  int i;
  int SelectedIndex;
  %TYPE% *Ranges;
  %TYPE% Random%TYPE%;
  %TYPE% TotalSum;
  %TYPE% CumulativeSum;

  // Calculate sum of elements

  TotalSum = (%TYPE%)0.0;
  for ( i = 0; i < NumberOfElements; i++ )
    TotalSum += vector->data[i * vector->stride];

  // Handle pathological case where sum of elements is zero

  if ( TotalSum == (%TYPE%)0.0 )
    {
      cout << "Vector%TYPE2%::CumulativeSumDraw() ERROR:";
      cout << "Cannot draw from vector where sum of elements is zero." << endl;
      exit( 1 );
    }

  // Ranges[] holds the cumulative sum per element so we can generate
  // a random number in the range 0.0..1.0 and find which bin it
  // falls into. Bins are bigger for larger elements in the array
   
  Ranges = new %TYPE%[NumberOfElements];

  // Calculate cumulative sum and store scaled version in Ranges[]

  CumulativeSum = (%TYPE%)0.0;
  for ( i = 0; i < NumberOfElements; i++ )
    {
      CumulativeSum += vector->data[i * vector->stride] / TotalSum;
      Ranges[i] = CumulativeSum;
    }

  // Random %TYPE% for bin selection

  Random%TYPE% = myrand();

  // SelectedIndex is the index corresponding to the Range[] bin that
  // brackets Random%TYPE%
   
  SelectedIndex = 0;
  for ( i = 1; i < NumberOfElements; i++ )
    if ( Random%TYPE% > Ranges[i-1] && Random%TYPE% <= Ranges[i] )
      SelectedIndex = i;

  delete [] Ranges;

  return( SelectedIndex );
}

void Vector%TYPE2%::Distinct()
{
  int flag;
  Vector%TYPE2% temp( NumberOfElements );

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

void Vector%TYPE2%::RemoveElement( int ElementNumber )
{
  Vector%TYPE2% temp;

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
    temp(i) = gsl_vector%GSL%_get( vector, i );

  gsl_vector%GSL%_free( vector );
  NumberOfElements--;
  vector = gsl_vector%GSL%_calloc( NumberOfElements );
  int counter = -1;
  for ( int i = 0; i < NumberOfElements + 1; i++ ){
    if( i != ElementNumber ){
      counter++;
      gsl_vector%GSL%_set( vector, counter, temp( i ) );
    }
  }   
}

void Vector%TYPE2%::AddElement( int ElementNumber )
{
  Vector%TYPE2% temp;

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
    temp(i) = gsl_vector%GSL%_get( vector, i );

  gsl_vector%GSL%_free( vector );
  NumberOfElements++;
  vector = gsl_vector%GSL%_calloc( NumberOfElements );
  int counter = 0;
  for ( int i = 0; i < NumberOfElements - 1; i++ ){
    if( i == ElementNumber ){
      gsl_vector%GSL%_set( vector, counter, 0 );
      counter++;
    }
    gsl_vector%GSL%_set( vector, counter, temp( i ) );
    counter++;
  }   
}

gsl_vector%GSL%* Vector%TYPE2%::GetGslVector()
{
  return( vector );
}
