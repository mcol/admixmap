#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include "rand.h"
#include "matrix_d.h"
#include "gaussian_d.h"

Gaussian::Gaussian()
{
   Dimension = 1;

   // Allocate vector containing mean

   Mean = new Matrix_d( Dimension, 1 );

   // Allocate matrix containing covariance matrix

   Covariance = new Matrix_d( Dimension, Dimension );
   InverseCovariance = new Matrix_d( Dimension, Dimension );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

Gaussian::Gaussian( int NewDimension )
{
   if ( NewDimension < 1 )
   {
      NewDimension = 1;
   }
   Dimension = NewDimension;

   // Allocate vector containing mean

   Mean = new Matrix_d( Dimension, 1 );

   // Allocate matrix containing covariance matrix

   Covariance = new Matrix_d( Dimension, Dimension );
   InverseCovariance = new Matrix_d( Dimension, Dimension );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

Gaussian::Gaussian( const Matrix_d &NewMean, const Matrix_d &NewCovariance )
{
   int NewDimension = NewMean.GetNumberOfRows();

   // Check dimensions of mean and covariance are equal
   
   if ( NewDimension != NewCovariance.GetNumberOfRows()
        || NewDimension != NewCovariance.GetNumberOfCols() )
   {
      std::cout << "ERROR: Unequal dimensions between mean and covariance." << std::endl;
      exit( 1 );
   }
   Dimension = NewDimension;

   // Allocate vector containing mean

   Mean = new Matrix_d( Dimension, 1 );
   *Mean = NewMean;

   // Allocate matrix containing covariance matrix

   Covariance = new Matrix_d( Dimension, Dimension );
   *Covariance = NewCovariance;
   InverseCovariance = new Matrix_d( Dimension, Dimension );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

Gaussian::Gaussian( const Gaussian &g )
{
   Dimension = g.GetDimension();

   // Allocate vector containing mean

   Mean = new Matrix_d( Dimension, 1 );

   // Allocate matrix containing covariance matrix

   Covariance = new Matrix_d( Dimension, Dimension );
   InverseCovariance = new Matrix_d( Dimension, Dimension );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;

   // Copy mean and covariance from g to *this

   *Mean = g.GetMean();
   *Covariance = g.GetCovariance();
}

Gaussian::~Gaussian()
{
   delete Mean;
   delete Covariance;
   delete InverseCovariance;
}

void Gaussian::SetDimension( int NewDimension )
{
   if ( Dimension == NewDimension )
      return;

   delete Mean;
   delete Covariance;
   delete InverseCovariance;

   if ( NewDimension < 1 )
      NewDimension = 1;
   Dimension = NewDimension;

   // Allocate vector containing mean

   Mean = new Matrix_d( Dimension, 1 );

   // Allocate matrix containing covariance matrix

   Covariance = new Matrix_d( Dimension, Dimension );
   InverseCovariance = new Matrix_d( Dimension, Dimension );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

Gaussian &Gaussian::operator=( const Gaussian &g )
{
   if ( &g == this )
      return *this;

   Dimension = g.GetDimension();

   delete Mean;
   delete Covariance;
   delete InverseCovariance;

   // Allocate vector containing mean

   Mean = new Matrix_d( Dimension, 1 );

   // Allocate matrix containing covariance matrix

   Covariance = new Matrix_d( Dimension, Dimension );
   InverseCovariance = new Matrix_d( Dimension, Dimension );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;

   // Copy mean and covariance from g to *this

   *Mean = g.GetMean();
   *Covariance = g.GetCovariance();

   return *this;
}

float Gaussian::GetCovariance( int Row, int Col ) const
{
   return( (*Covariance)( Row, Col ) );
}

Matrix_d &Gaussian::GetCovariance() const
{
   return( (*Covariance) );
}

Matrix_d &Gaussian::GetMean() const
{
   return( (*Mean) );
}

double &Gaussian::GetMean( int Element ) const
{
   return( (*Mean)( Element, 0 ) );
}

void Gaussian::SetMean( const Matrix_d &NewMean )
{
   *Mean = NewMean;
}

void Gaussian::SetCovariance( int row, int col, float NewValue )
{
   (*Covariance)(row,col) = NewValue;
   if ( row != col )
      (*Covariance)(col,row) = NewValue;
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

void Gaussian::SetCovariance( const Matrix_d &NewCovariance )
{
   *Covariance = NewCovariance;
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

void Gaussian::SetCovariance( double NewValue )
{
   Covariance->SetElements( NewValue );
   InverseCovarianceIsDirty = 1;
   CovarianceDeterminantIsDirty = 1;
}

std::ostream& operator<< ( std::ostream& os, const Gaussian& G )
{
   os << G.GetMean();
   os << G.GetCovariance();
//    os << G.GetMeanPrior();
//    os << G.GetCovariancePrior();

   return os;
}	

int Gaussian::GetDimension() const
{
   return( Dimension );
}

Matrix_d Gaussian::Draw()
{
   int i;
   Matrix_d draw( Dimension, 1 );
   Matrix_d cd;

   for ( i = 0; i < Dimension; i++ )
      draw( i, 0 ) = gennor( (double)0.0, (double)1.0 );

   if ( Covariance->CholeskyDecomposition( &cd ) == 0 )
   {
      std::cout << "Cholesky decomposition failed..." << std::endl;
      exit( 1 );
   }
   draw = ( cd * draw );
   draw += *Mean;
   
   return( draw );
}
