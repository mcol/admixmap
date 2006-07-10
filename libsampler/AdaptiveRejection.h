// *-*-C++-*-*
/** 
 *   AdaptiveRejection.h
 *   This class implements an adaptative rejection algorithm for  
 *   log-concave distributions. 
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ADAPTIVEREJECTION_H
#define ADAPTIVEREJECTION_H 1

#include <vector>
#include <math.h>
#include <stdlib.h>

/**
   object to hold the various attributes of a point used in ARS.
   Intended to be used in a vector and thus sortable
*/
class ARSPoint{
public:
  double abscissa;///<x
  double height;  ///< value of logdensity at x
  double gradient;  ///< derivative of logdensity at x
  double upper;  ///< value of upper hull at x
  double lower;  ///< value of lower hull at x
  double s;
  /**
     intersection of tangent at x[i] with tangent at x[i+1].
     last element of vector will have UpperBound as z
  */
  double z;
  ///area above x, ie under curve at tangent at x between z[i-1] and z[i], where z[-1] = LowerBound  
  double area;
  double cumarea;  //<cumulative area

  ///default constructor
  ARSPoint(){
    abscissa = 0.0; height = 0.0; gradient = 0.0; upper = 0.0; lower = 0.0; s = 0.0; z = 0.0; area = 0.0; cumarea = 0.0;
  }
///function for StrictWeakOrdering of points
  bool operator<(const ARSPoint& p) const {
    return (this->abscissa < p.abscissa);
  }

};

/**
 *   This class implements an adaptative rejection algorithm for  
 *   log-concave distributions. 
*/  
class AdaptiveRejection{

public:
  AdaptiveRejection();
  double Sample(const void* const args, double (*secondDeriv)(double, const void* const));
  double Sample(double mode, const void* const args, double (*secondDeriv)(double, const void* const));
  double ARS(const void* const args, double initialPoints[], unsigned numberOfInitialPoints);
  void Initialise(bool upperBound, bool lowerBound, double upper, double lower, 
		  double (*height) (double, const void* const),
		  double (*gradient)(double, const void* const));
  void setLowerBound(double);
  void setUpperBound(double); 

private:
  unsigned K;  ///<initial number of points
  std::vector< ARSPoint> Points;
  ARSPoint NewPoint;
  unsigned pos;  ///<position of new point
  bool hasUpperBound;
  bool hasLowerBound;
  double UpperBound;
  double LowerBound;
  double heightAtMode;

  double (*height) (double, const void* const);
  double (*gradient)(double, const void* const); 

  double AreaUnderTangentCurve(double x1, double h, double g, 
			       double z1, double z2, 
			       bool lower, bool upper); 
  double TransformPoint(double u, double g, double s1, double s2, double z1, double z2, 
			      bool lower, bool upper);

  double TangentIntersection(double x0, double x1, double h0, double h1, double g0, double g1);

  double UpperHull(double x, double x0, double h, double g);
  double LowerHull(double x, double x0, double x1, double h0, double h1, bool lower, bool upper);

  void InitialisePoints(double x[3], const void* const args);
  void SamplePoint(const void* const);
  bool TestNewPoint();
  void Update();
  void TestForLogConcavity();

  double SimpleModeSearch( const double a, const double b, const void* const args, 
			   double (*gradient)(double, const void* const) );

  double NewtonRaphson(const void* const args, double (*gradient)(double, const void* const), 
			       double (*secondDeriv)(double, const void* const) );

  void SetInitialPoints(double mode, double x[3], const void* const args, 
			double (*secondDeriv)(double, const void* const));

};

#endif











