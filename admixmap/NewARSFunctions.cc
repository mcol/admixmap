/** 
 *   ADMIXMAP
 *   DARS.cc
 *   This class implements an adaptative rejection algorithm for  
 *   log-concave distributions. 
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include "NewARS.h"
#include <iostream>
#include <fstream>

#define DEBUG 0 //set to 1 for debug info

using namespace std;

#define EPS  0.00000001            /* critical relative x-value difference */

NewDARS::NewDARS(){
  K = 3;//initial number of points
  hasUpperBound = false;
  hasLowerBound = false;
  UpperBound = 1.0;//irrelevent when no upper bound
  LowerBound = 0.0;//irrelevant when no lower bound
}

void NewDARS::Initialise(bool upperBound, bool lowerBound, double upper, double lower, 
			double (*h) (double, const double* const),
			double (*g)(double, const double* const)){
  hasUpperBound = upperBound;
  hasLowerBound = lowerBound;
  UpperBound = upper;
  LowerBound = lower;
  height = h;
  gradient = g;
}

//use this function if mode not known
double NewDARS::Sample(const double* const args, double (*secondDeriv)(double, const double* const)){
  K = 3;//number of initial points
  double x[] = {-2.0, 0.0, 2.0};//dummy args

  //code to find mode and initial points
  double dfa = 1, dfb = -1; //gradient at upper and lower bound (if there are any) should be +ve and -ve respectively
  
  if( hasLowerBound ) // if bounded below, evaluate gradient at lower bound 
    dfa = (*gradient)( LowerBound, args );
  if( hasUpperBound ) // if bounded above, evaluate gradient at upper bound
    dfb = (*gradient)( UpperBound, args );
  
  //mode not at boundary
  if( dfa > 0.0 && dfb < 0.0 ){
    // Mode search
    if( hasLowerBound && hasUpperBound ){ // if bounded below and above 
      SimpleModeSearch( LowerBound, UpperBound, x, args, gradient, secondDeriv );
    }
    else{ // if unbounded below or above
      NewtonRaphson(x, args, gradient, secondDeriv); // assigns x[1] as mode, x[0], x[2] as +/- 3 times 2nd deriv
    }
    
    // if log density infinite at x[0], assign x[0] as mean of x[0] and x[1]
    while( isinf( (*height)( x[0], args ) ) )
      x[0] = ( x[0] + x[1] ) / 2.0;
    // x[2] similarly
    while( isinf( (*height)( x[2], args ) ) )
      x[2] = ( x[2] + x[1] ) / 2.0;
  }

  else if( dfb > 0 ){ // mode at upper bound; 
    //only true if there is a UB since dfb initialise to -1
    x[2] = UpperBound;
    if( hasLowerBound ){ // if bounded below
      x[0] = LowerBound;
      x[1] = (x[0] + x[2]) / 2; // assign x[1] as mean of x[0] and x[2]
    }
    else{ // use gradient at max to assign x[0] and x[2]
      x[1] = x[2] - 2/dfb;
      x[0] = x[1] - 2/dfb;
    }
  }
  else{ // mode at lower bound
    //only true if there is an LB since dfa initialised to 1 
    x[0] = LowerBound;
    if( hasUpperBound ){ // assign x[1] as mean of x[0] and x[2]
      x[2] = UpperBound;
      x[1] = (x[0] + x[2]) / 2;
    }
    else{ // use gradient at min
      x[1] = x[0] - 2/dfa; 
      x[2] = x[1] - 2/dfa;
    }
  }

  return ARS(args, x, K);
}

//use this version if the mode is known but leaving object to determine two more points
// double NewDARS::ARS(const double* const args, double mode, double (*secondDeriv)(double, const double* const)){
//   K = 3;//number of initial points
//   //code for mode search goes here
//   //x[1] = .....

//   //use secondDeriv to determine x[0] and x[2]

// }

//use this function if initial points are known (to save time computing mode etc)
double NewDARS::ARS(const double* const args, double initialPoints[], unsigned numberOfInitialPoints){
  K = numberOfInitialPoints;
  InitialisePoints(initialPoints, args);

  TestForLogConcavity();
  //?? no need to Initialise sampling function s
  //?? no need to initialise lower hull ??

  bool success = false;
  do{
    //sample new point 
    SamplePoint(args);
    success = TestNewPoint();  //returns true if accepted

    if(!success){
      //add new point if rejected
      Update();
    }
  }
  while (!success);

  Points.clear();//delete ready for next time
  return NewPoint.abscissa;
}

//end of public interface
//***************************************************************************************************************

double NewDARS::AreaUnderTangentCurve(double x1, double h, double g, double z1, double z2, 
				     bool lower, bool upper) {
  // returns unnormalized area under exponential curve formed by tangent to log density  
  // arguments: x1 abscissa in the interval, h height of log density at x1, g gradient at x1, 
  // z1 and z2 lower and upper bounds ofinterval, 
  // lower=true if lower bound, upper=true if upper bound
  double unscaled = 0.0;
  // should check that z1 <= x1 if lower, and x1 <= z2 if upper
  // to prevent computational underflow, could return unscaled with h = log scale factor
  
  if(lower && upper){
    if(fabs(g) > EPS) {
      unscaled = (exp(g*(z2 - x1)) - exp(g*(z1 - x1))) / g;
    } else { // gradient close to 0
      // evaluate difference between exponentials using first two terms of series 
      unscaled = z2 - z1 + 0.5*(z2 + z1 - 2.0*x1)*(z2 - z1)*g;
    }
  } else if(lower && g < 0) {//no upper limit and negative gradient
    unscaled = -exp(g*(z1 - x1)) / g; 
  } else if(upper && g > 0) {//no lower limit and positive gradient
    unscaled = exp(g*(z2 - x1)) / g;
  } else {
    cout << "Error in arguments passed to AreaUnderTangentCurve \n";
    exit(1);
  }
  return exp(h) * unscaled;
} 

//turns a uniform point to a point on the x axis
double NewDARS::TransformPoint(double u, double x1, double g, double s1, double s2, double z1, double z2, 
			      bool lower, bool upper) {
  if(s1>s2 || u<s1 || u>s2 || (!lower && g<0) || (!upper && g>0) ){
    cerr<<"Error in arguments to TransformPoint"<<endl;
    exit(1);
  }
  double x;
  if(lower && upper){
    x = z1 + (z2-z1)*(u-s1) / (s2-s1);
  } else if(lower && g < 0) {//no upper limit and negative gradient
    x = z1 + (log(u-s1) - log(1.0 - s1)) / g;
  } else if(upper && g > 0) {//no lower limit and positive gradient
    x = z2 + (log(u) - log(s2)) / g;
  } else {
    cout << "Error in arguments passed to TransformPoint \n";
    exit(1);
  }
  return x;
} 

 
void NewDARS::InitialisePoints(double x[3], const double* const args){
  for(unsigned i = 0; i < K; ++i){
    ARSPoint p;
    p.abscissa = x[i];
    p.height = (*height)(x[i], args);
    p.gradient = (*gradient)(x[i], args);
    p.upper = p.height;//upper hull at x is the same as height at x
    Points.push_back(p);
  }
  stable_sort(Points.begin(), Points.end()); // PROBLEM: garbles points
  for(unsigned i = 0; i < K-1; ++i){
    Points[i].z = TangentIntersection(x[i], x[i+1], Points[i].height, Points[i+1].height, Points[i].gradient, Points[i+1].gradient);
  }
  Points[K-1].z = UpperBound;


}
  
void NewDARS::SamplePoint(const double* const args){
  //sample from standard uniform
  double u = myrand();
#if DEBUG ==1
  cout<<"Sampled u = "<<u<<endl;;
#endif
  
  //compute areas under tangents and hence interval probabilities
  //possibly wasteful to redo all each time
  double sum = 0.0;
  bool upper;
  Points[0].area = AreaUnderTangentCurve(Points[0].abscissa, Points[0].height, Points[0].gradient,
				  LowerBound, Points[0].z, hasLowerBound, true);
  sum += Points[0].area;

  for(unsigned i = 1; i < K; ++i){
    if(i == (K-1))upper = hasUpperBound; else upper=true;
    Points[i].area = AreaUnderTangentCurve(Points[i].abscissa, Points[i].height, Points[i].gradient, 
				    Points[i-1].z, Points[i].z, true, upper);
    sum += Points[i].area;
  }
  pos = 0;
  for(unsigned i = 0; i < K; ++i){
    Points[i].cumarea = 0.0;
    Points[i].area /= sum; //normalize areas
    for(unsigned j = 0; j <= i; ++j) Points[i].cumarea += Points[j].area;
    if(u > Points[i].cumarea)++pos;
  }
  //new point lies between z[pos-1] and z[pos] ie same interval as x[pos]

#if DEBUG ==1
  cout<<"pos = "<<pos<<endl;
  cout<<"x\tboundary\tCumArea"<<endl;
  for(unsigned i = 0; i < K; ++i){
    cout<<Points[i].abscissa<<"\t"<<Points[i].z<<"\t"<<Points[i].cumarea<<endl;
  }
#endif
  
  //transform from u to xnew
  bool lower = pos > 0 || hasLowerBound;
  upper = pos < K-1 || hasUpperBound;

  double z1, z2;
  double s1, s2;
  if(pos == 0) s1 = 0.0;
  else s1 = Points[pos-1].cumarea;
  if(pos==K-1)s2 = 1.0;
  else s2 = Points[pos].cumarea;
  if(pos>0)z1 = Points[pos-1].z; else z1 = LowerBound;
  z2 = Points[pos].z;
  
  NewPoint.abscissa = TransformPoint(u, Points[pos].abscissa, Points[pos].gradient, s1, s2, z1, z2, lower, upper);
  //check new point is in range
  if( ( (hasUpperBound || pos<K-1) && NewPoint.abscissa > Points[pos].z ) //too big
      || (pos>0 && NewPoint.abscissa < Points[pos-1].z) || (pos==0 && hasLowerBound && NewPoint.abscissa < LowerBound))//too small
    {
      cerr<<"Error: miscalculated point: "<<NewPoint.abscissa<<endl;
      exit(1);
    }
#if DEBUG ==1
  cout<<"New point at x = "<<NewPoint.abscissa<<endl;
#endif

  //compute values of functions at new point
  NewPoint.height = height(NewPoint.abscissa, args);//call to log density function
  NewPoint.gradient = gradient(NewPoint.abscissa, args);//call to gradient function
  NewPoint.upper = UpperHull(NewPoint.abscissa, Points[pos].abscissa, Points[pos].height, Points[pos].gradient);

  if(NewPoint.abscissa > Points[pos].abscissa) ++pos;//so that xnew < x[pos] 

  //compute lower hull
  lower =true; if(pos==0 && !hasLowerBound) lower=false;
  upper = true; if(pos==K && !hasUpperBound)upper=false;
  double h0 = 0.0;
  if(pos == 0 && hasLowerBound )h0 = (*height)(LowerBound, args);
  else if(pos > 0)h0 = Points[pos-1].height;
  double h1 = 0.0;
  if(pos==K && hasUpperBound)h1 = (*height)(UpperBound, args);
  else if(pos < K)h1 = Points[pos].height;

 
  NewPoint.lower = LowerHull(NewPoint.abscissa, Points[pos-1].abscissa, Points[pos].abscissa, 
			     Points[pos-1].height, Points[pos].height, lower, upper);
}

bool NewDARS::TestNewPoint(){
  double w = myrand();
  bool accept = false;
  //squeeze test
  if (w <= exp(NewPoint.lower - NewPoint.upper) ) accept = true;
  else
    //rejection test
    if(w <= exp(NewPoint.height - NewPoint.upper) ) accept = true;
  
  return accept;
}

void NewDARS::Update(){
  //height and gradient at new point already calculated
  //as are upper and lower hulls
  Points.insert(Points.begin()+pos, 1, NewPoint);//insert before pos
#if DEBUG ==1
  cout<<"new point added"<<endl<<endl;
#endif

  //increment dimension
  ++K;

  //re-sort the Points
  //stable_sort(Points.begin(), Points.end());

  TestForLogConcavity();

  //recalculate intersection points. Note that not all need be recalculated
  unsigned i0 = 0;
  if(pos > 0)i0 = pos-1;
  for(unsigned i = i0; i < K-1; ++i){
    Points[i].z = TangentIntersection(Points[i].abscissa, Points[i+1].abscissa, Points[i].height, Points[i+1].height, 
				      Points[i].gradient, Points[i+1].gradient);
  }
  Points[K-1].z = UpperBound;
}

double NewDARS::TangentIntersection(double x0, double x1, double h0, double h1, double g0, double g1){
  //computes x coord of point of intersection of tangents at x0 and x1
  //h0, h1 = heights at two points
  //g0, g1 = gradients at the two points
  if(fabs(g0 - g1) < EPS){
    cerr<<"Error in args passed to TangentIntersection: g0 = "<<g0<<", g1 = "<<g1<<endl;
    exit(1);
  }
  return (h1 - h0 - x1*g1 + x0*g0) / (g0 - g1);
}

double NewDARS::UpperHull(double x, double x0, double h, double g){
  //computes value of upper hull at x
  //x0 = the point in the same interval as x
  //h, g = height and gradient at x0
  return h + (x-x0)*g;
}
double NewDARS::LowerHull(double x, double x0, double x1, double h0, double h1, bool lower, bool upper){
  //computes value of lower hull between x0 and x1 at x
  if( (lower && x < x0) || (upper && x>x1)){
    cerr<<"Error in LowerHull: "<<x0<<" "<<x<<" "<<x1<<endl;exit(1);
  }
  double lh;
  if(upper && lower){
    lh = (h0*(x1-x) + h1*(x-x0)) / (x1-x0);
  }
  else if(upper){//no lower bound
    lh = h1;
  }
  else if(lower){//no upper bound
    lh = h0;
  }
  return lh;
}

void NewDARS::TestForLogConcavity(){
  for(unsigned i = 0; i < K-1; ++i){
    if(Points[i].gradient < Points[i+1].gradient){
      cerr<<"Error: failed test for log-concavity"<<endl;
      exit(2);
    }
  }
}

void NewDARS::SimpleModeSearch( double aa, double bb, double *x, const double* const args,
			       double (*gradient)(double, const double* const), double (*secondDeriv)(double, const double* const) )
//searches for mode between aa and bb
//places mode and a point either side in x
{
  double newnum;
  double a, b, num, df2, df1, dfa, dfb, ddf;
  int count = 0;
  a = aa + 0.00000001;
  b = bb - 0.00000001;
  
  newnum = ( a + b ) / 2;
  do{
    count++;
    df1 = (*gradient)(newnum, args);
    if( df1 > 0.0 ){
      num = (newnum + b) / 2;
      df2 = (*gradient)(num, args);
      if( df2 < df1 ){
	a = newnum;
	newnum = num;
	df1 = df2;
      }
      else{
	b = num;
      }
    }
    else{
      num = (newnum + a) / 2;
      df2 = (*gradient)(num, args);
      if( df2 > df1 ){
	b = newnum;
	newnum = num;
	df1 = df2;
      }
      else{
	a = num;
      }
    }
  }while( fabs(df1) > 0.01 );
  dfa = (*gradient)(a, args);
  dfb = (*gradient)(b, args);
  x[1] = newnum;
  
  if( dfa > 1.0 )
    x[0] = a;
  else{
    ddf = (*secondDeriv)( newnum, args );
    x[0] = x[1] + 3.0 / ddf;
    if( hasLowerBound && x[0] < LowerBound )
      x[0] = LowerBound;
  }
  
  if( dfb < -1.0 )
    x[2] = b;
  else{
    ddf = (*secondDeriv)( newnum, args );
    x[2] = x[1] - 3.0 / ddf;
    if( hasUpperBound && x[2] > UpperBound )
      x[2] = UpperBound;
  }

}

void NewDARS::NewtonRaphson(double *x, const double* const args, double (*gradient)(double, const double* const), 
			     double (*secondDeriv)(double, const double* const) )
{
// Using modified Newton-Raphson (step size * 2) to find two points
// either side of the mode. Given these two points use most simple and
// robust mode search SimpleModeSearch(). Two points chosen such that
// they are 3 * second derivative from the mode.
  double newnum, oldnum, step, dfnew, dfold, ddf, a, b;
  
  do{
    oldnum = newnum;
    ddf = (*secondDeriv)(oldnum, args );
    dfold = (*gradient)(oldnum, args );
    step = -dfold / ddf;
    newnum += 2 * step;
    if( hasLowerBound && newnum < LowerBound )
      newnum = LowerBound;
    if( hasUpperBound && newnum > UpperBound )
      newnum = UpperBound;
    dfnew = (*gradient)(newnum, args );
    // continues looping until oldnum and new num are on opposite sides of mode, or gradient zero at newnum  
  }while( ( (dfold * dfnew) > 0 ) && fabs(dfnew) > 0.01 ); 
  
  if( fabs(dfnew) > 0.01 ){ // if gradient not zero at newnum, set a to lower value and b to higher value
    if( oldnum < newnum ){
      a = oldnum;
      b = newnum;
    }
    else{
      a = newnum;
      b = oldnum;
    }
    // if a below lower bound or b below upper bound, reset to bounds 
    if( hasLowerBound && a < LowerBound )
      a = LowerBound;
    else if( hasUpperBound && b > UpperBound )
      b = UpperBound;
    
    SimpleModeSearch( a, b, x, args, gradient, secondDeriv );
  }
  else{
    x[1] = newnum;
    x[0] = x[1] + 3.0 / ddf;
    x[2] = x[1] - 3.0 / ddf;
    if( hasLowerBound && x[0] < LowerBound )
      x[0] = LowerBound;
    else if( hasUpperBound && x[2] > UpperBound )
      x[2] = UpperBound;
  }
}

// ********* test code ***************

void NewDARS::test(){
  // test algorithm
  double x1 = 3.0;
  double h =  -10.0;
  double z1 = -5.0;
  double z2 = 15.0;
  double gradient = 1.5;
  bool lower = false;
  bool upper=true;

  //cout << AreaUnderTangentCurve(x1, h, gradient, z1, z2, lower, upper) << "\n";
  cout<<AreaUnderTangentCurve(-1, -0.5, 1.0, 0.0, -0.5, false, true)<<endl;
  cout<<AreaUnderTangentCurve(0, 0, 0, -0.5, 0.5, true, true)<<endl;
  cout<<AreaUnderTangentCurve(1, -0.5, -1.0, 0.5, 0.0, true, false)<<endl;
}

double LogNormalDensity(double x, const double* const args){
  double mean = args[0];
  double sd = args[1];

  return -0.5 * log(sd) - 0.5 * (x - mean)*(x-mean) / sd;
}
double LogNormalGradient(double x, const double* const args){
  double mean = args[0];
  double sd = args[1];

  return (mean - x) / sd;
}
double LogNormalSecondDeriv(double x, const double* const args){
  double sd = args[1];
  return -1.0 / sd;
}


int main() {
  NewDARS Sampler;
  ofstream outfile("./ARSPoints.txt");//outputfile to inspect points

  //test algorithm with univariate Normal distribution
  Sampler.Initialise(false, false, 0, 0, LogNormalDensity, LogNormalGradient);

  double args[2] = {0.0, 1.0};// mean, sd
  for(unsigned i = 0; i < 1000; ++i){
    double x  = Sampler.Sample(args, LogNormalSecondDeriv);
    cout<<"Sampled point "<<i+1<< ": "<<x<<endl<<endl;
    outfile<<x<<" ";
  }

  outfile.close();
  //Sampler.test(); 
}
