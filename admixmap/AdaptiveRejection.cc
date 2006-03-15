/** 
 *   ADMIXMAP
 *   AdaptiveRejection.cc
 *   This class implements an adaptative rejection algorithm for  
 *   log-concave distributions. 
 *   Copyright (c) 2005 LSHTM
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
#include "AdaptiveRejection.h"
#include "rand.h"
#include <iostream>


#define DEBUG 0 //set to 1 for debug info

using namespace std;

#define EPS  0.00000001            /* critical relative x-value difference */

AdaptiveRejection::AdaptiveRejection(){
  K = 3;//initial number of points
  hasUpperBound = false;
  hasLowerBound = false;
  UpperBound = 1.0;//irrelevent when no upper bound
  LowerBound = 0.0;//irrelevant when no lower bound
}

void AdaptiveRejection::Initialise(bool upperBound, bool lowerBound, double upper, double lower, 
			double (*h) (double, const void* const),
			double (*g)(double, const void* const)){
  hasUpperBound = upperBound;
  hasLowerBound = lowerBound;
  UpperBound = upper;
  LowerBound = lower;
  height = h;
  gradient = g;
}

void AdaptiveRejection::setLowerBound(double LB){
  LowerBound = LB;
  hasLowerBound = true;
}
void AdaptiveRejection::setUpperBound(double UB){
  UpperBound = UB;
  hasUpperBound = true;
}

//use this function if mode not known
double AdaptiveRejection::Sample(const void* const args, double (*secondDeriv)(double, const void* const)){
  K = 3;//number of initial points
  double x[] = {-2.0, 0.0, 2.0};//dummy args

  //code to find mode and initial points
  double dfa = 1, dfb = -1; //gradient at upper and lower bound (if there are any) should be +ve and -ve respectively
  double mode;  

  if( hasLowerBound ) // if bounded below, evaluate gradient just above lower bound 
    dfa = (*gradient)( LowerBound +EPS, args );
  //if(isinf(dfa))dfa = (*gradient)(LowerBound+EPS, args);
  if( hasUpperBound ) // if bounded above, evaluate gradient just below upper bound
    dfb = (*gradient)( UpperBound - EPS, args );
  //if(isinf(dfb))dfb = (*gradient)(UpperBound-EPS, args);
  
  //mode not at boundary
  if( dfa > 0.0 && dfb < 0.0 ){
    // Mode search
    if( hasLowerBound && hasUpperBound ){ // if bounded below and above 
      mode = SimpleModeSearch(LowerBound, UpperBound, args, gradient);
    }
    else{ // if unbounded below or above
      mode = NewtonRaphson(args, gradient, secondDeriv); // assigns x[1] as mode, x[0], x[2] as +/- 3 times 2nd deriv
    }
    SetInitialPoints(mode, x, args, secondDeriv);    

    // if log density infinite at x[0], assign x[0] as mean of x[0] and x[1]
    while( isinf( (*height)( x[0], args ) ) )
      x[0] = ( x[0] + x[1] ) / 2.0;
    // x[2] similarly
    while( isinf( (*height)( x[2], args ) ) )
      x[2] = ( x[2] + x[1] ) / 2.0;
  }

  else if( dfb > 0 ){ // mode at upper bound; 
    //only true if there is a UB since dfb initialised to -1
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
    //only true if there is a LB since dfa initialised to 1 
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
 double AdaptiveRejection::Sample(double mode, const void* const args, double (*secondDeriv)(double, const void* const)){
   K = 3;//number of initial points
  double x[] = {-2.0, 0.0, 2.0};//dummy args

  double dfa = 1, dfb = -1; //gradient at upper and lower bound (if there are any) should be +ve and -ve respectively
  
  if( hasLowerBound ) // if bounded below, evaluate gradient at lower bound 
    dfa = (*gradient)( LowerBound, args );
  //if(isinf(dfa))dfa = (*gradient)(LowerBound+EPS, args);
  if( hasUpperBound ) // if bounded above, evaluate gradient at upper bound
    dfb = (*gradient)( UpperBound, args );
  //if(isinf(dfb))dfb = (*gradient)(UpperBound-EPS, args);
  
  //mode not at boundary
  if( dfa > 0.0 && dfb < 0.0 ){
    SetInitialPoints(mode, x, args, secondDeriv);    

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

//use this function if initial points are known (to save time computing mode etc)
double AdaptiveRejection::ARS(const void* const args, double initialPoints[], unsigned numberOfInitialPoints){
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

double AdaptiveRejection::AreaUnderTangentCurve(double x1, double h, double g, double z1, double z2, 
				     bool lower, bool upper) {
  // returns unnormalized area under exponential curve formed by tangent to log density  
  // arguments: x1 abscissa in the interval, h height of log density at x1, g gradient at x1, 
  // z1 and z2 lower and upper bounds ofinterval, 
  // lower=true if lower bound, upper=true if upper bound
  double unscaled = 0.0;
  if(z2 < z1 ){
    cerr<<"AdaptiveRejectionSampler: Error in arguments passed to AreaUnderTangentCurve: z1 > z2 \n"; 
    exit(1);
  }
  if(lower && (x1 < z1)){
    cerr<<"AdaptiveRejectionSampler: Error in arguments passed to AreaUnderTangentCurve: x < z1 \n"; 
    exit(1);
  }
  if(upper && (x1 > z2)){
    cerr<<"AdaptiveRejectionSampler: Error in arguments passed to AreaUnderTangentCurve: x > z2 \n"; 
    exit(1);
  }

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
    cerr << "AdaptiveRejectionSampler: Error in arguments passed to AreaUnderTangentCurve \n";
    exit(1);
  }
  return exp(h) * unscaled;
} 

//turns a uniform point to a point on the x axis
double AdaptiveRejection::TransformPoint(double u, double g, double s1, double s2, double z1, double z2, 
			      bool lower, bool upper) {
  //u = draw from standard uniform, g = gradient of tangent at point in interval, s1, s2 = cumulative areas at z1 and z2
  //z1, z2 = endpoints of interval on real line, lower, upper = indicators for lower and upper bounds
  if(s1>s2 || u<s1 || u>s2 || (!lower && g<0) || (!upper && g>0) ){
    cerr<<"Error in arguments to TransformPoint"<<endl;
    exit(1);
  }
  double x = u;
  double prop = (u-s1) / (s2-s1);//proportion of area between z1 and z2 that is to the left of u
  double maxterm = max(log(1.0 - prop) + g*z1, log(prop)+g*z2);//in case very steep gradient
  if(lower && upper){
    if(fabs(g) > EPS)
    x = (maxterm + log( exp( (g*z1) + log(1.0 - prop) - maxterm) + exp((g*z2) + log(prop)-maxterm) )) / g;
    else //if g close to zero, we effectively have to transform from U(s1,s2) to U(z1, z2)
      x = z1 + (z2 - z1) * prop;
  } else if(lower && g < 0) {//no upper limit and negative gradient
    x = z1 + exp( log(-log( 1.0 - prop)) - log( -g ) );
  } else if(upper && g > 0) {//no lower limit and positive gradient
    x = z2 - exp( log(-log(prop)) - log( g ) );
  } else {
    cout << "AdaptiveRejectionSampler: Error in arguments passed to TransformPoint \n";
    exit(1);
  }
  return x;
} 

 
void AdaptiveRejection::InitialisePoints(double x[3], const void* const args){
  sort(x, x+3);
  heightAtMode = 0.0;
  for(unsigned i = 0; i < K; ++i){
    ARSPoint p;
    p.abscissa = x[i];
    p.height = (*height)(x[i], args);
    //if(p.height < heightAtMode) heightAtMode = p.height;
    p.gradient = (*gradient)(x[i], args);
    p.upper = p.height;//upper hull at x is the same as height at x
    Points.push_back(p);
  }
  heightAtMode = Points[1].height;
  for(unsigned i = 0; i < K; ++i){
    Points[i].height -= heightAtMode;
    Points[i].upper -= heightAtMode;
  }

  stable_sort(Points.begin(), Points.end()); // PROBLEM: can garble points; redundant if x is sorted
  for(unsigned i = 0; i < K-1; ++i){
    Points[i].z = TangentIntersection(Points[i].abscissa, Points[i+1].abscissa, Points[i].height, Points[i+1].height, 
				      Points[i].gradient, Points[i+1].gradient);
  }
  Points[K-1].z = UpperBound;


}
  
void AdaptiveRejection::SamplePoint(const void* const args){
  //sample from standard uniform
  double u = Rand::myrand();
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
  cout<<"x         boundary  Area    Cum.Area  height  gradient"<<endl;
  for(unsigned i = 0; i < K; ++i){
    cout<<Points[i].abscissa<<"  "<<Points[i].z<<"  "<<Points[i].area<<"  "<<Points[i].cumarea<<"  "<<Points[i].height<<"  "<<Points[i].gradient<<endl;
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
  
  NewPoint.abscissa = TransformPoint(u, Points[pos].gradient, s1, s2, z1, z2, lower, upper);
  //check new point is in range
  if( ( (hasUpperBound || pos<K-1) && NewPoint.abscissa > Points[pos].z ) //too big
      || (pos>0 && NewPoint.abscissa < Points[pos-1].z) || (pos==0 && hasLowerBound && NewPoint.abscissa < LowerBound))//too small
    {
      cerr<<"AdaptiveRejectionSampler: Error: miscalculated point: "<<NewPoint.abscissa<<endl;
      exit(1);
    }
#if DEBUG ==1
  cout<<"New point at x = "<<NewPoint.abscissa<<endl;
#endif

  //compute values of functions at new point
  NewPoint.height = height(NewPoint.abscissa, args);
  NewPoint.height -= heightAtMode;
  NewPoint.gradient = gradient(NewPoint.abscissa, args);
  NewPoint.upper = UpperHull(NewPoint.abscissa, Points[pos].abscissa, Points[pos].height, Points[pos].gradient);

  if(NewPoint.abscissa > Points[pos].abscissa) ++pos;//so that xnew < x[pos] 

  //compute lower hull
  lower =true; if(pos==0 && !hasLowerBound) lower=false;
  upper = true; if(pos==K && !hasUpperBound)upper=false;
  double h0 = 0.0;
  if(pos == 0 && hasLowerBound )h0 = (*height)(LowerBound+EPS, args);
  else if(pos > 0)h0 = Points[pos-1].height;
  double h1 = 0.0;
  if(pos==K && hasUpperBound)h1 = (*height)(UpperBound-EPS, args);
  else if(pos < K)h1 = Points[pos].height;
  double x0 = LowerBound, x1=UpperBound;
  if(pos > 0)x0 = Points[pos-1].abscissa;
  if(pos < K)z1 = Points[pos].abscissa;
  NewPoint.lower = LowerHull(NewPoint.abscissa, x0, x1, h0, h1, lower, upper);
}

bool AdaptiveRejection::TestNewPoint(){
  double w = Rand::myrand();
  bool accept = false;
  //squeeze test
  if (w <= exp(NewPoint.lower - NewPoint.upper) ) {accept = true;
  }
  else{
    //rejection test
    if(w <= exp(NewPoint.height - NewPoint.upper) ) accept = true;
    }
  
  return accept;
}

void AdaptiveRejection::Update(){

  //insert new point in array before pos
  Points.insert(Points.begin()+pos, 1, NewPoint);
  //increment dimension
  ++K;
 
#if DEBUG ==1
  cout<<"new point added"<<endl<<endl;
#endif

  //re-sort the Points
  //stable_sort(Points.begin(), Points.end());

  TestForLogConcavity();

  //height, gradient, upper and lower hulls at new point already calculated

  //recalculate intersection points. Note that not all need be recalculated
  unsigned i0 = 0;
  if(pos > 0)i0 = pos-1;
  for(unsigned i = i0; i < K-1; ++i){
    Points[i].z = TangentIntersection(Points[i].abscissa, Points[i+1].abscissa, Points[i].height, Points[i+1].height, 
				      Points[i].gradient, Points[i+1].gradient);
  }
  Points[K-1].z = UpperBound;
}

double AdaptiveRejection::TangentIntersection(double x0, double x1, double h0, double h1, double g0, double g1){
  //computes x coord of point of intersection of tangents at x0 and x1
  //h0, h1 = heights at two points
  //g0, g1 = gradients at the two points
  if(fabs(g0 - g1) < EPS){
    cerr<<"AdaptiveRejectionSampler: Error in args passed to TangentIntersection: g0 = "<<g0<<", g1 = "<<g1<<endl;
    exit(1);
  }
  return (h1 - h0 - x1*g1 + x0*g0) / (g0 - g1);
}

double AdaptiveRejection::UpperHull(double x, double x0, double h, double g){
  //computes value of upper hull at x
  //x0 = the point in the same interval as x
  //h, g = height and gradient at x0
  return h + (x-x0)*g;
}
double AdaptiveRejection::LowerHull(double x, double x0, double x1, double h0, double h1, bool lower, bool upper){
  //computes value of lower hull between x0 and x1 at x
  if( (lower && x < x0) || (upper && x>x1)){
    cerr<<"Error in LowerHull: "<<x0<<" "<<x<<" "<<x1<<endl;exit(1);
  }
  double lh = 0.0;
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

void AdaptiveRejection::TestForLogConcavity(){
  for(unsigned i = 0; i < K-1; ++i){
    if(Points[i].gradient < Points[i+1].gradient){
      cerr<<"AdaptiveRejectionSampler: Error: failed test for log-concavity"<<endl;
      for(unsigned j = 0; j < K; ++j)
	cerr<<"x["<<j<<"] = "<<Points[j].abscissa<<", g["<<j<<"] = "<<Points[j].gradient<<endl;
      exit(2);
    }
  }
}

double AdaptiveRejection::SimpleModeSearch( const double a, const double b, const void* const args,
					    double (*gradient)(double, const void* const) )
  //searches for mode between a and b, without evaluating gradient at a or b 
{
  double x = 0.0, gradx = 0.0, upper = 0.0, lower = 0.0, gradupper = 0.0, gradlower = 0.0, 
    deriv2nd = 0.0;
  int count = 0;
  // step 1: find two values either side of mode at which gradient can be evaluated
  x = 0.5 * (a + b);
  gradx = (*gradient)(x, args);
  if(gradx > 0.0) { // step to right until gradient is negative 
    lower = x;
    gradlower = gradx;
    do {
      x = 0.5 * (x + b);
      gradx = (*gradient)(x, args);
      ++count;
    } while(gradx > 0.0);
    upper = x;
    gradupper = gradx;
  } else {
    upper = x;
    gradupper = gradx;
    do {
      x = 0.5 * (a + x);
      gradx = (*gradient)(x, args);
      ++count;
    } while(gradx < 0.0);
    lower = x;
    gradlower = gradx;
  }
  // step 2: iterate using average 2nd derivative for approximate newton-raphson update 
  do {
    deriv2nd = (gradupper - gradlower) / (upper - lower);
    x = lower - gradlower / deriv2nd;
    gradx = (*gradient)(x, args);
    if(gradx > 0.0) {
      lower = x;
      gradlower = gradx;
    } else {
      upper = x;
      gradupper = gradx;
    }
    ++count;
  } while(fabs(gradx) > 0.001);
  return x;
}

double AdaptiveRejection::NewtonRaphson(const void* const args, double (*gradient)(double, const void* const),
			    double (*secondDeriv)(double, const void* const) )
{
// Using modified Newton-Raphson (step size * 2) to find two points
// either side of the mode. Given these two points use most simple and
// robust mode search SimpleModeSearch(). Two points chosen such that
// they are 3 * second derivative from the mode.
  double newnum = 0.0, oldnum, step, dfnew, dfold, ddf, a, b;
  
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
  
  double mode;
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
    
    mode = SimpleModeSearch( a, b, args, gradient);

  }//end if(fabs(dfnew > 0.01)
  else{
    mode = newnum;
  }
  return mode;
}

void AdaptiveRejection::SetInitialPoints(double mode, double x[3], const void* const args, 
			       double (*secondDeriv)(double, const void* const)){
  //given a mode and function for 2nd derivative, determines two points either side of the mode and places them in x

  x[1] = mode;
  
  double ddf = (*secondDeriv)( mode, args );
  x[0] = x[1] - 2.5 / sqrt(-ddf);
  if( hasLowerBound && x[0] < LowerBound )
    x[0] = LowerBound;

  ddf = (*secondDeriv)( mode, args );
  x[2] = x[1] + 2.5 / sqrt(-ddf);
  if( hasUpperBound && x[2] > UpperBound )
    x[2] = UpperBound;
}
