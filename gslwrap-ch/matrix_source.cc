#include "rand.h"

#include "vector.h"
#include "vector_d.h"
#include "vector_i.h"

#include "matrix.h"
#include "matrix_d.h"
#include "matrix_i.h"

#define MATRIX_TYPE%TYPE2%
#define BOUNDS_CHECK 1

using namespace std;

Matrix%TYPE2%::Matrix%TYPE2%():
  _matrix((gsl_matrix%GSL%*)malloc(sizeof(gsl_matrix%GSL%)))
{
  assert(_matrix);
  _matrix->size1 = 0;
  _matrix->size2 = 0;
  _matrix->tda   = 0;
  _matrix->data  = 0;
  _matrix->block = 0;
  _matrix->owner = 0;
}

Matrix%TYPE2%::Matrix%TYPE2%(size_t NewRows, size_t NewCols):
  _matrix(gsl_matrix%GSL%_calloc(NewRows,NewCols))
{
  assert(NewRows != 0);
  assert(NewCols != 0);
  assert(_matrix);
}

Matrix%TYPE2%::Matrix%TYPE2%(const Matrix%TYPE2%& other):
  _matrix(gsl_matrix%GSL%_calloc(other.GetNumberOfRows(),other.GetNumberOfCols()))
{
  _MissingRow = other._MissingRow;
  _MissingCol = other._MissingCol;
  gsl_matrix%GSL%_memcpy(_matrix,other._matrix);
}

Matrix%TYPE2%::~Matrix%TYPE2%()
{
  gsl_matrix%GSL%_free( _matrix );
  _matrix = 0;
}

Matrix%TYPE2% &Matrix%TYPE2%::operator=( const Matrix& m )
{
  if((void*)this == (void*)&m){
    return *this;
  }

  if (GetNumberOfRows() != m.GetNumberOfRows() ||
      GetNumberOfCols() != m.GetNumberOfCols()) {
    gsl_matrix%GSL%_free( _matrix );
    _matrix = gsl_matrix%GSL%_calloc(m.GetNumberOfRows(),m.GetNumberOfCols());
  }
  
  _MissingRow.clear();
  _MissingCol.clear();
  for (int i = 0; i < m.GetNumberOfMissingValues(); i++ ){
    _MissingRow.push_back(m.GetMissingValueRow(i));
    _MissingCol.push_back(m.GetMissingValueCol(i));
  }

  for(int i=0; i<GetNumberOfRows(); i++){
    for(int j=0; j<GetNumberOfCols(); j++){
      (*this)(i,j) = static_cast<%TYPE%>(m(i,j));
    }
  }
  
  return *this;
}

Matrix%TYPE2% &Matrix%TYPE2%::operator=( const Matrix_d& m )
{
  if((void*)this == (void*)&m){
    return *this;
  }

  if (GetNumberOfRows() != m.GetNumberOfRows() ||
      GetNumberOfCols() != m.GetNumberOfCols()) {
    gsl_matrix%GSL%_free( _matrix );
    _matrix = gsl_matrix%GSL%_calloc(m.GetNumberOfRows(),m.GetNumberOfCols());
  }
  
  _MissingRow.clear();
  _MissingCol.clear();
  for (int i = 0; i < m.GetNumberOfMissingValues(); i++ ){
    _MissingRow.push_back(m.GetMissingValueRow(i));
    _MissingCol.push_back(m.GetMissingValueCol(i));
  }

  for(int i=0; i<GetNumberOfRows(); i++){
    for(int j=0; j<GetNumberOfCols(); j++){
      (*this)(i,j) = static_cast<%TYPE%>(m(i,j));
    }
  }
  
  return *this;
}

Matrix%TYPE2% &Matrix%TYPE2%::operator=( const Matrix_i& m )
{
  if((void*)this == (void*)&m){
    return *this;
  }

  if (GetNumberOfRows() != m.GetNumberOfRows() ||
      GetNumberOfCols() != m.GetNumberOfCols()) {
    gsl_matrix%GSL%_free( _matrix );
    _matrix = gsl_matrix%GSL%_calloc(m.GetNumberOfRows(),m.GetNumberOfCols());
  }
  
  _MissingRow.clear();
  _MissingCol.clear();
  for (int i = 0; i < m.GetNumberOfMissingValues(); i++ ){
    _MissingRow.push_back(m.GetMissingValueRow(i));
    _MissingCol.push_back(m.GetMissingValueCol(i));
  }

  for(int i=0; i<GetNumberOfRows(); i++){
    for(int j=0; j<GetNumberOfCols(); j++){
      (*this)(i,j) = static_cast<%TYPE%>(m(i,j));
    }
  }
  
  return *this;
}

%TYPE%&
Matrix%TYPE2%::operator()( int row, int col )
{
  return _matrix->data[row * GetNumberOfCols() + col];
}

%TYPE%&
Matrix%TYPE2%::operator()( int row, int col ) const
{
  return _matrix->data[row * GetNumberOfCols() + col];
}

Matrix_d Matrix%TYPE2%::Double() const
{
  Matrix_d m(GetNumberOfRows(),GetNumberOfCols());

  for(int i=0; i<GetNumberOfRows(); i++){
    for(int j=0; j<GetNumberOfCols(); j++){
      m(i,j) = (double)(*this)(i,j);
    }
  }
  return(m);
}

Matrix Matrix%TYPE2%::Float() const
{
  Matrix m(GetNumberOfRows(),GetNumberOfCols());

  for(int i=0; i<GetNumberOfRows(); i++){
    for(int j=0; j<GetNumberOfCols(); j++){
      m(i,j) = (float)(*this)(i,j);
    }
  }
  return(m);
}

Matrix_i Matrix%TYPE2%::Integer() const
{
  Matrix_i m(GetNumberOfRows(),GetNumberOfCols());

  for(int i=0; i<GetNumberOfRows(); i++){
    for(int j=0; j<GetNumberOfCols(); j++){
      m(i,j) = (int)(*this)(i,j);
    }
  }
  return(m);
}

int Matrix%TYPE2%::GetNumberOfRows() const
{
  return _matrix->size1;
}

int Matrix%TYPE2%::GetNumberOfCols() const
{
  return _matrix->size2;
}

void Matrix%TYPE2%::SetNumberOfElements(size_t NewRows, size_t NewCols)
{
  _MissingRow.clear();
  _MissingCol.clear();

  assert( NewRows != 0 );
  assert( NewCols != 0 );
  if(NewRows==0 || NewCols==0){
    cerr << endl << endl << endl;
    cerr << "BUG - Matrix%TYPE2%::SetNumberOfElements(size_t,size_t) was passed zero as an argument... but continuing" << endl;
    cerr << "BUG -   first arg==" << NewRows << "; second arg==" << NewCols << endl;
    if(NewRows==0){
      NewRows=1;
      cerr << "FIXME: setting first arg to 1" << endl;
    }
    if(NewCols==0){
      NewCols=1;
      cerr << "FIXME: setting second arg to 1" << endl;
    }
    cerr << endl << endl << endl;
  }


  if (GetNumberOfRows() != (int)NewRows ||
      GetNumberOfCols() != (int)NewCols ){
    gsl_matrix%GSL%_free( _matrix );
    _matrix = gsl_matrix%GSL%_calloc( NewRows, NewCols );
  } else {
    SetElements(static_cast<%TYPE%>(0.0));
  }
}

void Matrix%TYPE2%::SetElements(%TYPE% NewValue)
{
  for (int i=0; i<GetNumberOfRows(); i++){
    for (int j=0; j<GetNumberOfCols(); j++){
      (*this)(i,j) = NewValue;
    }
  }
}

void Matrix%TYPE2%::SetDiagonal(%TYPE% NewValue)
{
#ifdef BOUNDS_CHECK
  assert(GetNumberOfRows()==GetNumberOfCols());
#endif
  for (int i=0; i<GetNumberOfRows(); i++){
    (*this)(i,i) = NewValue;
  }
}

void Matrix%TYPE2%::Load(const char *FileName)
{
  FileName >> *this;
}

ostream& operator<< (ostream& os,const Matrix%TYPE2%& M)
{
  int i;
  int j;

  os.setf(ios::fixed);

  for (i=0; i<M.GetNumberOfRows(); i++){
    for (j=0; j<M.GetNumberOfCols(); j++){
      if (M.IsMissingValue(i,j)){
	os << setw( 11 ) << "#" << " ";
      } else {
	os << setprecision( 6 ) << setw( 11 ) << M(i,j) << " ";
      }
    }
    os << endl;
  }
  return os;
}	

void operator>> ( const char * FileName, Matrix%TYPE2%& M )
{
  int NumRows = 0;
  int NumCols = 0;
  int NumMissingValues = 0;
  int PrevNumCols = 0;
  int MissingElementsThisRow = 0;
  int CurrentMissingValue = 0;
  ifstream fileOne;
  ifstream fileTwo;
  Vector%TYPE2% CurrentRow;

  M._MissingRow.clear();
  M._MissingCol.clear();

  fileOne.open(FileName);
  assert(fileOne.is_open());
  while( !fileOne.eof() ){
    fileOne >> CurrentRow;
    
    NumMissingValues += CurrentRow.GetNumberOfMissingValues();
    NumCols = CurrentRow.GetNumberOfElements();
    
#ifdef BOUNDS_CHECK
    // Check for ragged input
    assert(NumRows==0 || NumCols==PrevNumCols);
#endif
    PrevNumCols = NumCols;
    
    NumRows++;
  }
  fileOne.close();
  
  M.SetNumberOfElements(NumRows,NumCols);
  M._MissingRow.assign(NumMissingValues,0);
  M._MissingCol.assign(NumMissingValues,0);

  // Rewind to beginning of file and
  // reset eof flag
   
   
  fileTwo.open(FileName);
   
  // Load matrix elements

  CurrentMissingValue = 0;
  for (int i = 0; i < M.GetNumberOfRows(); i++ ){
    fileTwo >> CurrentRow;
    M.SetRow( i, CurrentRow );
    MissingElementsThisRow = CurrentRow.GetNumberOfMissingValues();
    
    // If there are missing elements store their row and column indices
    
    if ( MissingElementsThisRow ){
      for (int j = 0; j < CurrentRow.GetNumberOfElements(); j++ ){
	if ( CurrentRow.IsMissingElement( j ) ){
	  M._MissingRow[CurrentMissingValue] = i;
	  M._MissingCol[CurrentMissingValue] = j;
	  CurrentMissingValue++;
	}
      }
    }
  }
  
  fileTwo.close();
  return;
}	

double Matrix%TYPE2%::Determinant() const
{
  int sign;
  double LogDeterminant;

  LogDeterminant = this->LogDeterminant( &sign );

  return( sign * exp( LogDeterminant ) );
}

double Matrix%TYPE2%::LogDeterminant( int *sign ) const
{
  // Perform lower-upper decomposition on a copy of this matrix
  Matrix%TYPE2% temp = (*this);
  gsl_permutation *permutation = gsl_permutation_calloc(GetNumberOfRows());
  *sign = temp.LU_decomp(permutation);
  gsl_permutation_free( permutation );

  // return sum of logarithms of the diagonal elements
  return temp.LU_lndet();
}

void Matrix%TYPE2%::InvertUsingLUDecomposition()
{
  int ReturnValue;
  int sign;
  gsl_matrix *LowerUpper;
  gsl_matrix *Inverse;
  gsl_permutation *permutation;

  // Copy matrix into lu because gsl_la_decomp_LU_impl destroys the
  // matrix it is sent. GSL checks whether matrix is square etc.
  // Can't use gsl_matrix_copy because matrix is possibly double and
  // gsl_matrix copy assumes that it's double.
  Inverse = gsl_matrix_calloc(GetNumberOfRows(),GetNumberOfCols());
  LowerUpper = gsl_matrix_calloc(GetNumberOfRows(),GetNumberOfCols());
  permutation = gsl_permutation_calloc(GetNumberOfRows());
  for (int i=0; i<GetNumberOfRows(); i++){
    for (int j=0; j<GetNumberOfCols(); j++){
      int offset = i * Inverse->size2 + j;
      Inverse->data[offset] = _matrix->data[offset];
      LowerUpper->data[offset] = _matrix->data[offset];
    }
  }

  // Perform lower-upper decomposition
  ReturnValue = gsl_linalg_LU_decomp( LowerUpper, permutation, &sign );
  assert(ReturnValue==GSL_SUCCESS);

  // Calculate inverse using LU decomposed form of matrix
  gsl_linalg_LU_invert( LowerUpper, permutation, Inverse );

  // Copy results from Inverse back into our %TYPE% matrix[][]
  for (int i=0; i<GetNumberOfRows(); i++){
    for (int j=0; j<GetNumberOfCols(); j++){
      int offset = i * Inverse->size2 + j;
      _matrix->data[offset] = static_cast<%TYPE%>(Inverse->data[offset]);
    }
  }

  // Free memory
  gsl_permutation_free( permutation );
  gsl_matrix_free( LowerUpper );
  gsl_matrix_free( Inverse );
}

Vector%TYPE2% Matrix%TYPE2%::GetRow( int row ) const
{
  Vector%TYPE2% temp(GetNumberOfCols());
  for(int col=0; col<GetNumberOfCols(); col++){
    temp(col) = (*this)(row,col);
  }
  return temp;
}

Vector%TYPE2% Matrix%TYPE2%::GetColumn( int col ) const
{
  Vector%TYPE2% temp(GetNumberOfRows());
  for(int row=0; row<GetNumberOfRows(); row++){
    temp(row) = (*this)(row,col);
  }
  return( temp );
}

void Matrix%TYPE2%::SetRow(int row,const Vector%TYPE2%& v)
{
#ifdef BOUNDS_CHECK
  assert(GetNumberOfCols() == v.GetNumberOfElements());
  assert(row<GetNumberOfRows());
#endif   
  for (int i=0;i<GetNumberOfCols();i++){
    (*this)(row,i) = v(i);
  }
}

void Matrix%TYPE2%::SetColumn(int col,const Vector%TYPE2%& v)
{
#ifdef BOUNDS_CHECK
  assert(GetNumberOfRows() == v.GetNumberOfElements());
  assert(col<GetNumberOfCols());
#endif   
  for (int i=0;i<GetNumberOfRows();i++){
    (*this)(i,col) = v(i);
  }
}

void Matrix%TYPE2%::Symmetrize()
{
  %TYPE% m;

#ifdef BOUNDS_CHECK
  assert(GetNumberOfRows() == GetNumberOfCols());
#endif   
  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = i + 1; j < GetNumberOfCols(); j++ )
	{
	  m = (*this)(i,j) + (*this)(j,i);
	  m /= 2;
	  (*this)(i,j) = (*this)(j,i) = m;
	}
    }
}

Vector%TYPE2% Matrix%TYPE2%::GetDiagonal() const
{
  Vector%TYPE2% temp( GetNumberOfRows() );

#ifdef BOUNDS_CHECK
  assert(GetNumberOfRows() == GetNumberOfCols());
#endif   

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      temp(i) = (*this)(i,i);
    }
  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::SubMatrix( int rowstart, int rowfinish, int colstart, int colfinish ) const
{
#ifdef BOUNDS_CHECK
  assert(rowstart < GetNumberOfRows());
  assert(rowfinish < GetNumberOfRows());
  assert(rowstart >= 0);
  assert(rowfinish >= 0);
  assert(colstart < GetNumberOfCols());
  assert(colfinish < GetNumberOfCols());
  assert(colstart >= 0);
  assert(colfinish >= 0);
  assert(rowstart <= rowfinish);
  assert(colstart <= colfinish);
#endif

  Matrix%TYPE2% temp( rowfinish - rowstart + 1, colfinish - colstart + 1 );

  for ( int row = rowstart; row <= rowfinish; row++ ){
    for ( int col = colstart; col <= colfinish; col++ ){
      temp( row - rowstart, col - colstart ) = (*this)(row,col);
    }
  }
  return( temp );
}

void Matrix%TYPE2%::SubMatrix2( int rowstart, int rowfinish, int colstart, int colfinish )
{
#ifdef BOUNDS_CHECK
  assert(rowstart < GetNumberOfRows());
  assert(rowfinish < GetNumberOfRows());
  assert(rowstart >= 0);
  assert(rowfinish >= 0);
  assert(colstart < GetNumberOfCols());
  assert(colfinish < GetNumberOfCols());
  assert(colstart >= 0);
  assert(colfinish >= 0);
  assert(rowstart <= rowfinish);
  assert(colstart <= colfinish);
#endif
  Matrix%TYPE2% temp;
  temp.SetNumberOfElements(rowfinish - rowstart + 1, colfinish - colstart + 1 );
  int count = 0, count2 = 0;
  Vector%TYPE2% rowmissing( GetNumberOfMissingValues() ), colmissing( GetNumberOfMissingValues() );

  for ( int row = rowstart; row <= rowfinish; row++ ){
    for ( int col = colstart; col <= colfinish; col++ ){
      temp( row - rowstart, col - colstart ) = (*this)(row,col);
      if( count2 < GetNumberOfMissingValues() ){
	while( (int)_MissingRow[ count2 ] < row || 
	       ( (int)_MissingRow[ count2 ] == row && (int)_MissingCol[ count2 ] < col ) )
	  count2++;
	if( (int)_MissingRow[ count2 ] == row && (int)_MissingCol[ count2 ] == col ){
	  rowmissing( count ) = row - rowstart;
	  colmissing( count ) = col - colstart;
	  count++;
	  count2++;
	}
      }
    }
  }

  (*this) = temp;

  _MissingRow.assign(count,0);
  _MissingCol.assign(count,0);
  // Copy the missing value row and column indices if necessary
  for ( int i = 0; i < GetNumberOfMissingValues(); i++ ){
    (int)_MissingRow[i] = rowmissing(i);
    (int)_MissingCol[i] = colmissing(i);
  }

}

Matrix%TYPE2% ConcatenateHorizontally( const Matrix%TYPE2% &m1, const Matrix%TYPE2% &m2 )
{
#ifdef BOUNDS_CHECK
  assert(m1.GetNumberOfRows()==m2.GetNumberOfRows());
#endif

  Matrix%TYPE2% temp( m1.GetNumberOfRows(), m1.GetNumberOfCols() + m2.GetNumberOfCols() );

  for (int col = 0; col < m1.GetNumberOfCols(); col++ )
    temp.SetColumn( col, m1.GetColumn( col ) );
  for (int col = 0; col < m2.GetNumberOfCols(); col++ )
    temp.SetColumn( col + m1.GetNumberOfCols(), m2.GetColumn( col ) );

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator+( const Matrix &m ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  assert( GetNumberOfRows() == m.GetNumberOfRows() && GetNumberOfCols() == m.GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) + (%TYPE%)m(i,j);
	}
    }

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator+( const Matrix_d &m ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  assert( GetNumberOfRows() == m.GetNumberOfRows() && GetNumberOfCols() == m.GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) + (%TYPE%)m(i,j);
	}
    }

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator+( const Matrix_i &m ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  assert( GetNumberOfRows() == m.GetNumberOfRows() && GetNumberOfCols() == m.GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) + (%TYPE%)m(i,j);
	}
    }

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator+( %TYPE% f ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) + f;
	}
    }

  return( temp );
}

Matrix%TYPE2% &Matrix%TYPE2%::operator+=( const Matrix%TYPE2% &m )
{
  assert( GetNumberOfRows() == m.GetNumberOfRows() && GetNumberOfCols() == m.GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  (*this)(i,j) += m(i,j);
	}
    }

  return( *this );
}

Matrix%TYPE2% &Matrix%TYPE2%::operator*=( const %TYPE% f )
{
  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  (*this)(i,j) *= f;
	}
    }

  return( *this );
}

Matrix%TYPE2% Matrix%TYPE2%::operator-( const Matrix%TYPE2% &m ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  assert( GetNumberOfRows() == m.GetNumberOfRows() && GetNumberOfCols() == m.GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) - m(i,j);
	}
    }

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator*( const Matrix%TYPE2% &m ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), m.GetNumberOfCols() );

#ifdef BOUNDS_CHECK
  assert(GetNumberOfCols()==m.GetNumberOfRows());
#endif

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < m.GetNumberOfCols(); j++ )
	{
	  %TYPE% d = (%TYPE%)0.0;
	  for (int k = 0; k < m.GetNumberOfRows(); k++ )
	    d += (*this)(i,k) * m(k,j);
	  temp(i,j) = d;
	}
    }

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator*( %TYPE% f ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) * f;
	}
    }

  return( temp );
}

Matrix%TYPE2% operator*( %TYPE% f, const Matrix%TYPE2% &m )
{
  Matrix%TYPE2% temp( m.GetNumberOfRows(), m.GetNumberOfCols() );

  for (int i = 0; i < m.GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < m.GetNumberOfCols(); j++ )
	{
	  temp(i,j) = m(i,j) * f;
	}
    }

  return( temp );
}

Matrix%TYPE2% Matrix%TYPE2%::operator/( %TYPE% f ) const
{
  Matrix%TYPE2% temp( GetNumberOfRows(), GetNumberOfCols() );

  assert(f != 0.0);
  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(i,j) = (*this)(i,j) / f;
	}
    }

  return( temp );
}

Matrix%TYPE2% &Matrix%TYPE2%::operator/=( %TYPE% f )
{
  assert(f != 0.0);
  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  (*this)(i,j) /= f;
	}
    }

  return( *this );
}

int Matrix%TYPE2%::GetMissingValueRow( int index ) const
{
#ifdef BOUNDS_CHECK
  assert(index < GetNumberOfMissingValues());
  assert(index >= 0);
#endif
  return( _MissingRow[index] );
}

int Matrix%TYPE2%::GetMissingValueCol( int index ) const
{
#ifdef BOUNDS_CHECK
  assert(index < GetNumberOfMissingValues());
  assert(index >= 0);
#endif
  return( _MissingCol[index] );
}

int Matrix%TYPE2%::GetNumberOfMissingValues() const
{
  return _MissingRow.size();
}

bool Matrix%TYPE2%::IsMissingValue( int row, int col ) const
{
#ifdef BOUNDS_CHECK
  assert(row < GetNumberOfRows());
  assert(row >= 0);
  assert(col < GetNumberOfCols());
  assert(col >= 0);
#endif
  for (int i = 0; i < GetNumberOfMissingValues(); i++ ){
    if ( row == (int)_MissingRow[i] && col == (int)_MissingCol[i] ){
      return true;
    }
  }
  return false;
}

Vector%TYPE2% Matrix%TYPE2%::ColumnSumExcludingMissingValues() const
{
  Vector%TYPE2% temp( GetNumberOfCols() );

  temp = ColumnSum();
  for (int i = 0; i < GetNumberOfMissingValues(); i++ ){
    temp( _MissingCol[i] ) -= (*this)(_MissingRow[i],_MissingCol[i]);
  }
  return( temp );
}

Vector%TYPE2% Matrix%TYPE2%::ColumnMeanExcludingMissingValues() const
{
  Vector%TYPE2% NonMissingValueCount( GetNumberOfCols() );

  // Count non-missing values for each column

  NonMissingValueCount.SetElements( (%TYPE%)GetNumberOfRows() );
  for (int i = 0; i < GetNumberOfMissingValues(); i++ )
    {
      NonMissingValueCount( _MissingCol[i] )--;
    }
  // Check whether all values in a column are missing to avoid
  // division by zero

#ifdef BOUNDS_CHECK
  for (int i = 0; i < GetNumberOfCols(); i++ ){
    assert(NonMissingValueCount( i ) != 0.0);
  }
#endif
  return( ColumnSumExcludingMissingValues() / NonMissingValueCount );
}

void Matrix%TYPE2%::SetMissingValuesToColumnMeans()
{
  Vector%TYPE2% MissingValueEstimate = ColumnMeanExcludingMissingValues();

  for (int i = 0; i < GetNumberOfMissingValues(); i++ ){
    (*this)(_MissingRow[i],_MissingCol[i]) = MissingValueEstimate( _MissingCol[i] );
  }
}

Matrix%TYPE2% Matrix%TYPE2%::Transpose() const
{
  Matrix%TYPE2% temp( GetNumberOfCols(), GetNumberOfRows() );

  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  temp(j,i) = (*this)(i,j);
	}
    }

  return( temp );
}

int Matrix%TYPE2%::CholeskyDecomposition( Matrix%TYPE2% *CD ) const
{
  %TYPE% sum;
  %TYPE% **a;
  Vector%TYPE2% p;

  assert(GetNumberOfRows() == GetNumberOfCols());
  assert(IsSymmetric());

  // Construct Cholesky factor matrix

  CD->SetNumberOfElements(GetNumberOfRows(),GetNumberOfRows());

  // Allocate temporary %TYPE% working matrix and copy contents of matrix[][]

  a = new %TYPE% * [GetNumberOfRows()];
  for (int i = 0; i < GetNumberOfRows(); i++ )
    a[i] = new %TYPE%[GetNumberOfCols()];
  for (int i = 0; i < GetNumberOfRows(); i++ )
    for (int j = 0; j < GetNumberOfCols(); j++ )
      a[i][j] = (*this)(i,j);

  // Construct Vector%TYPE2% for Cholesky diagonal

  p.SetNumberOfElements( GetNumberOfRows() );

  int k;
  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = i; j < GetNumberOfRows(); j++ )
	{
	  for ( sum = a[i][j], k = ((int)i)-1; k >= 0; k-- )
            sum -= a[i][k] * a[j][k];
	  if ( i == j )
	    {
	      if ( sum <= 0.0 )
		{
		  if ( fabs((double)sum ) < 1.0e-6 )
		    sum *= (%TYPE%)-1.0;
		  else
		    {
		      for ( i = 0; i < GetNumberOfRows(); i++ )
			delete a[i];
		      delete [] a;
		      return( 0 );
		    }
		}
	      p( i ) = (%TYPE%) sqrt((double)sum );
	    }
	  else
            if ( p( i ) != 0.0 )
	      a[j][i] = sum / p( i );
	}
    }

  for (int i = 0; i < GetNumberOfRows(); i++ ){
    for (int j = 0; j < GetNumberOfCols(); j++ ){
      if (i==j){
	(*CD)(i,i) = p(i);
      } else if (i>j) {
	(*CD)(i,j) = a[i][j];
      } else {
	(*CD)(i,j) = static_cast<%TYPE%>(0.0);
      }
    }
  }
  // Free working matrix a[][]
   
  for (int i = 0; i < GetNumberOfRows(); i++ )
    delete a[i];
  delete [] a;

  return( 1 );
}

Vector%TYPE2% Matrix%TYPE2%::ColumnSum() const
{
  Vector%TYPE2% temp( GetNumberOfCols() );

  temp.SetElements( (%TYPE%)0.0 );
  for (int i = 0; i < GetNumberOfRows(); i++ )
    for (int j = 0; j < GetNumberOfCols(); j++ )
      temp(j) += (*this)(i,j);

  return( temp );
}

Vector%TYPE2% Matrix%TYPE2%::RowSum() const
{
  Vector%TYPE2% temp( GetNumberOfRows() );

  temp.SetElements( (%TYPE%)0.0 );
  for (int i = 0; i < GetNumberOfCols(); i++ )
    for (int j = 0; j < GetNumberOfRows(); j++ )
      temp(j) += (*this)(j,i);

  return( temp );
}

Vector%TYPE2% Matrix%TYPE2%::ColumnMean() const
{
  return( ColumnSum() / GetNumberOfRows() );
}

Vector%TYPE2% Matrix%TYPE2%::RowMean() const
{
  return( RowSum() / GetNumberOfCols() );
}

bool Matrix%TYPE2%::IsSymmetric() const
{
  for (int i = 0; i < GetNumberOfRows(); i++ ){
    for (int j = i+1; j < GetNumberOfCols(); j++ ){
      if ((*this)(i,j) != (*this)(j,i)){
	return false;
      }
    }
  }
  return true;
}

void Matrix%TYPE2%::Eigenvalue2( Vector%TYPE2% *Eigenvalue )
{
  gsl_matrix *GSLData;
  gsl_vector *GSLVal;
  gsl_eigen_symm_workspace *w;

  GSLData = gsl_matrix_calloc( GetNumberOfRows(), GetNumberOfCols() );
  for (int i = 0; i < GetNumberOfRows(); i++ )
    {
      for (int j = 0; j < GetNumberOfCols(); j++ )
	{
	  int offset = i * GSLData->size2 + j;
	  GSLData->data[offset] = _matrix->data[offset];
	}
    }
  GSLVal = gsl_vector_calloc( GetNumberOfRows() );
  w = gsl_eigen_symm_alloc( GetNumberOfRows() );
  
  assert(gsl_eigen_symm( GSLData, GSLVal, w ) == GSL_SUCCESS);
  
  Eigenvalue->SetNumberOfElements( GetNumberOfRows() );
  for(int i = 0; i < GetNumberOfRows(); i++ )
    {
      (*Eigenvalue)(i) = ((%TYPE%)(GSLVal->data[i * GSLVal->stride]));
    }

  gsl_eigen_symm_free( w ); 
  gsl_vector_free( GSLVal );
  gsl_matrix_free( GSLData );

  return;
}

// This function finds the indexes of the elements
// of Vector%TYPE2% whose values are equal to t
// Returns a two-column matrix, where the first column
// stores the rows, the second the columns.
Matrix%TYPE2% Matrix%TYPE2%::FindElement(%TYPE% t)
{
  int n = 0;
  Matrix%TYPE2% Index;
  Matrix%TYPE2% Tmp(GetNumberOfRows()*GetNumberOfCols(), 2);

  Tmp.SetElements((%TYPE%)0.0);
   
  for (int i = 0; i < GetNumberOfRows(); i++ )
    for (int j = 0; j < GetNumberOfCols(); j++ ){
      if ( (*this)(i,j) == t ){
	Tmp(n,0) = i;
	Tmp(n,1) = j;
	n++;
      }
    }

  Index.SetNumberOfElements(n,2);

  for (int i = 0; i < n; i++ ){
    Index(i,0) = Tmp(i,0);
    Index(i,1) = Tmp(i,1);
  }

  return(Index);
}

/**
 * Computes the logarithm of the absolute value of the determinant of 
 * a matrix from its LU decomposition (*this). This function may be
 * useful is the direct computation of the determinant would underflow
 * or overflow.
 */
double
Matrix%TYPE2%::LU_lndet()
{
  double lndet = 0.0;
  for (int i=0; i<GetNumberOfRows(); i++){
    lndet += log (fabs(static_cast<double>((*this)(i,i))));
  }
  return lndet;
}

/**
 * Factorizes the square matrix into the LU decomposition, given a 
 * permutation matrix.
 */
int
Matrix%TYPE2%::LU_decomp (gsl_permutation * p){
  assert(GetNumberOfRows() == GetNumberOfCols());
  assert((int)p->size == GetNumberOfRows());

  int signum = 1;
  gsl_permutation_init (p);

  for (int j=0; j < GetNumberOfRows()-1; j++){
    /* Find maximum in the j-th column */

    %TYPE% ajj, max = static_cast<%TYPE%>(fabs(static_cast<double>((*this)(j,j))));
    size_t i_pivot = j;

    for (int i = j + 1; i < GetNumberOfRows(); i++){
      %TYPE% aij = static_cast<%TYPE%>(fabs(static_cast<double>((*this)(i,j))));

      if (aij > max){
	max = aij;
	i_pivot = i;
      }
    }

    if ((int)i_pivot != j){
      gsl_matrix%GSL%_swap_rows (_matrix, j, i_pivot);
      gsl_permutation_swap (p, j, i_pivot);
      signum = -signum;
    }

    ajj = (*this)(j,j);

    if (ajj != 0.0){
      for (int i = j + 1; i < GetNumberOfRows(); i++){
	%TYPE% aij = (*this)(i,j) / ajj;
	(*this)(i,j) = aij;

	for (int k = j + 1; k < GetNumberOfRows(); k++){
	  %TYPE% aik = (*this)(i,k);
	  %TYPE% ajk = (*this)(j,k);
	  (*this)(i,k) = aik - aij * ajk;
	}
      }
    }
  }
      
  return signum;
}

void Matrix%TYPE2%::SetMissingElement(int row, int col)
{
    if (!IsMissingValue(row, col)) {
        _MissingRow.push_back(row);
        _MissingCol.push_back(col);
    }
}
