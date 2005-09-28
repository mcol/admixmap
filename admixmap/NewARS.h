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
  
double myrand()

/* to return a standard uniform random number */
{
   return ((double)rand() + 0.5)/((double)RAND_MAX + 1.0);
}

class NewARS{

public:
  NewARS();
  double ARS(const double* const);
  void Initialise(bool upperBound, bool lowerBound, double upper, double lower, 
		  double (*height) (double, const double* const),
		  double (*gradient)(double, const double* const)); 

  void test();//temporary, just for demo purposes

private:
  unsigned K;//initial number of points
  std::vector< ARSPoint> Points;
  ARSPoint NewPoint;
  bool hasUpperBound;
  bool hasLowerBound;
  double UpperBound;
  double LowerBound;

  double (*height) (double, const double* const);
  double (*gradient)(double, const double* const); 

  double AreaUnderTangentCurve(double x1, double h, double g, 
			       double z1, double z2, 
			       bool lower, bool upper); 
  double TransformPoint(double u, double x1, double g, double s, double z1, double z2, 
			bool lower, bool upper);

  double TangentIntersection(double x0, double x1, double h0, double h1, double g0, double g1);

  double UpperHull(double x, double x0, double h, double g);

  void InitialisePoints(double x[3], const double* const args);
  void SamplePoint(const double* const);
  bool Sample(const double* const args);
  void Update();


};

#endif











