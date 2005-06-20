// gaussian.h

#if !defined( _gaussian_h )

#define _gaussian_h

#include <fstream>

///
class Matrix_d;

///
class Gaussian
{
public:
   ///
   Gaussian();
   ///
   Gaussian( int NewDimension );
   ///
   Gaussian( const Gaussian &g );
   ///
   Gaussian( const Matrix_d &NewMean, const Matrix_d &NewCovariance );
   ///
   ~Gaussian();
   ///
   void SetDimension( int NewDimension );
   ///
   Gaussian &operator=( const Gaussian &g );
   ///
   void Randomize();
   ///
   float GetCovariance( int Row, int Col ) const;
   ///
   Matrix_d &GetMean() const;
   ///
   double &GetMean( int Element ) const;
   ///
   Matrix_d &GetCovariance() const;
   ///
   int GetDimension() const;
   ///
   friend std::ostream& operator<< ( std::ostream& os, const Gaussian& G );
   ///
   void SetMean( int index, float NewValue );
   ///
   void SetMean( const Matrix_d &NewMean );
   ///
   void AddToMean( const Matrix_d &a );
   ///
   void UpdateMean( const Matrix_d &a, int PriorCount );
   ///
   void SetCovariance( int row, int col, float NewValue );
   ///
   void SetCovariance( const Matrix_d &NewCovariance );
   ///
   void SetCovariance( double NewValue );
   ///
   double Likelihood( const Matrix_d &d );
   ///
   double LogLikelihood( const Matrix_d &d );
   ///
   Matrix_d Draw();
   void Draw(double *);
   ///
   ///
private:
   ///
   int Dimension;
   ///
   int InverseCovarianceIsDirty;
   ///
   int CovarianceDeterminantIsDirty;
   ///
   float CovarianceDeterminant;
   ///
   float NormalisationConstant;
   ///
   Matrix_d *Mean;
   ///
   Matrix_d *Covariance;
   ///
   Matrix_d *InverseCovariance;
   ///
};

#endif
