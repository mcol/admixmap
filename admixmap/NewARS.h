// *-*-C++-*-*
#ifndef NEWARS_H
#define NEWARS_H 1

#include <vector>
#include <math.h>
#include <stdlib.h>

//object to hold the various attributes of a point used in ARS
//intended to be used in a vector and thus sortable
class ARSPoint{
public:
  double abscissa;//x
  double height;//value of  logdensity at x
  double gradient;//derivative of logdensity at x
  double upper;//value of upper hull at x
  double lower;//value of lower hull at x
  double s;
  double z;//intersection of tangent at x_i with tangent at x_(i+1)
  // last element of vector will UpperBound as z
  double area;//area above x, ie under curve at tangent at x between z[i-1] and z[i], where z[-1] = LowerBound 
  double cumarea;//cumulative area

  //default constructor
  ARSPoint(){
    abscissa = 0.0; height = 0.0; gradient = 0.0; upper = 0.0; lower = 0.0; s = 0.0; z = 0.0; area = 0.0; cumarea = 0.0;
  }
//function for StrictWeakOrdering of points
  bool operator<(const ARSPoint& p) const {
    return (this->abscissa < p.abscissa);
  }

};
  
class NewDARS{

public:
  NewDARS();
  double Sample(const double* const args, double (*secondDeriv)(double, const double* const));
  double Sample(double mode, const double* const args, double (*secondDeriv)(double, const double* const));
  double ARS(const double* const args, double initialPoints[], unsigned numberOfInitialPoints);
  void Initialise(bool upperBound, bool lowerBound, double upper, double lower, 
		  double (*height) (double, const double* const),
		  double (*gradient)(double, const double* const));
  void setLowerBound(double);
  void setUpperBound(double); 

  void test();//temporary, just for demo purposes

private:
  unsigned K;//initial number of points
  std::vector< ARSPoint> Points;
  ARSPoint NewPoint;
  unsigned pos;//position of new point
  bool hasUpperBound;
  bool hasLowerBound;
  double UpperBound;
  double LowerBound;
  double heightAtMode;

  double (*height) (double, const double* const);
  double (*gradient)(double, const double* const); 

  double AreaUnderTangentCurve(double x1, double h, double g, 
			       double z1, double z2, 
			       bool lower, bool upper); 
  double TransformPoint(double u, double g, double s1, double s2, double z1, double z2, 
			      bool lower, bool upper);

  double TangentIntersection(double x0, double x1, double h0, double h1, double g0, double g1);

  double UpperHull(double x, double x0, double h, double g);
  double LowerHull(double x, double x0, double x1, double h0, double h1, bool lower, bool upper);

  void InitialisePoints(double x[3], const double* const args);
  void SamplePoint(const double* const);
  bool TestNewPoint();
  void Update();
  void TestForLogConcavity();

  double SimpleModeSearch( double aa, double bb, const double* const args, double (*gradient)(double, const double* const) );

  double NewtonRaphson(const double* const args, double (*gradient)(double, const double* const), 
			       double (*secondDeriv)(double, const double* const) );

  void SetInitialPoints(double mode, double x[3], const double* const args, 
			double (*secondDeriv)(double, const double* const));

};

#endif











