#include "NewARS.h"
#include <iostream>

using namespace std;

#define EPS  0.00000001            /* critical relative x-value difference */

NewARS::NewARS(){
  K = 3;//initial number of points
  hasUpperBound = false;
  hasLowerBound = false;
  UpperBound = 1.0;//irrelevent when no upper bound
  LowerBound = 0.0;//irrelevant when no lower bound
}

void NewARS::Initialise(bool upperBound, bool lowerBound, double upper, double lower, 
			double (*h) (double, const double* const),
			double (*g)(double, const double* const)){
  hasUpperBound = upperBound;
  hasLowerBound = lowerBound;
  UpperBound = upper;
  LowerBound = lower;
  height = h;
  gradient = g;
}

double NewARS::AreaUnderTangentCurve(double x1, double h, double g, double z1, double z2, 
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
double NewARS::TransformPoint(double u, double x1, double g, double s, double z1, double z2, 
			      bool lower, bool upper) {
  double x;
  if(lower && upper){
    if(fabs(g) > EPS) {
      x = ( log( u / (exp(z2*g) - exp(z1*g)) + exp(z1*g)   ) ) / g;
    } else { // gradient close to 0
      x = z1 + (u - s) / (z2-z1);//needs fixing, only valid for g==0
    }
  } else if(lower && g < 0) {//no upper limit and negative gradient
      x = z2 + log(u) / g;
  } else if(upper && g > 0) {//no lower limit and positive gradient
      x = z1 + log(u) / g;
  } else {
    cout << "Error in arguments passed to TransformPoint \n";
    exit(1);
  }
  return x;
} 

 
void NewARS::InitialisePoints(double x[3], const double* const args){
  for(unsigned i = 0; i < K; ++i){
    ARSPoint p;
    p.abscissa = x[i];
    p.height = (*height)(x[i], args);
    p.gradient = (*gradient)(x[i], args);
    p.upper = p.height;//upper hull at x is the same as height at x
    Points.push_back(p);
  }
  for(unsigned i = 0; i < K-1; ++i){
    Points[i].z = TangentIntersection(x[i], x[i+1], Points[i].height, Points[i+1].height, Points[i].gradient, Points[i+1].gradient);
  }
  Points[K-1].z = UpperBound;

  sort(Points.begin(), Points.end()); 
}
  
void NewARS::SamplePoint(const double* const args){
  //sample from standard uniform
  double u = 0.193;//myrand();
  cout<<"Sampled u = "<<u<<endl;;
  
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
  unsigned pos = 0;
  for(unsigned i = 0; i < K; ++i){
    Points[i].cumarea = 0.0;
    Points[i].area /= sum; //normalize areas
    for(unsigned j = 0; j <= i; ++j) Points[i].cumarea += Points[j].area;
    if(u > Points[i].cumarea)++pos;
  }
  //new point lies between z[pos-1] and z[pos] ie same interval as x[pos]
  
  //transform from x to xnew
  bool lower = pos > 0 || hasLowerBound;
  upper = pos < K-1 || hasUpperBound;
  double z1, z2;
  double s;
  if(pos == 0) s = 0.0;
  else s = Points[pos-1].cumarea;
  if(pos>0)z1 = Points[pos-1].z; else z1 = LowerBound;
  z2 = Points[pos].z;
  NewPoint.abscissa = TransformPoint(u, Points[pos].abscissa, Points[pos].gradient, s, z1, z2, lower, upper);

  cout<<"New point at x = "<<NewPoint.abscissa<<", pos= "<<pos<<endl;

  //compute values of functions at new point
  NewPoint.height = height(NewPoint.abscissa, args);//call to log density function
  NewPoint.gradient = gradient(NewPoint.abscissa, args);//call to gradient function
  NewPoint.upper = UpperHull(NewPoint.abscissa, Points[pos].abscissa, Points[pos].height, Points[pos].gradient);
  //NewPoint.lower = LowerHull(...);
  unsigned above = (NewPoint.abscissa < Points[pos+1].abscissa)? pos+1 : pos+2;//index of x above new point
  NewPoint.z = TangentIntersection(NewPoint.abscissa, Points[above].abscissa, NewPoint.height, Points[above].height,
				   NewPoint.gradient, Points[above].height);

}

bool NewARS::Sample(const double* const args){
  // sample new point 
  SamplePoint(args);

  double w = myrand();
  bool accept = false;
  //squeeze test
  if (w <= exp(NewPoint.lower - NewPoint.upper) ) accept = true;
  else
    //rejection test
    if(w <= exp(NewPoint.height - NewPoint.upper) ) accept = true;
  
  return accept;
}

void NewARS::Update(){
  //height and gradient at new point already calculated
  //as are upper and lower hulls
  Points.push_back(NewPoint);
  cout<<"new point added"<<endl;

  //increment dimension
  ++K;

  //re-sort the Points
  sort(Points.begin(), Points.end());
}

double NewARS::TangentIntersection(double x0, double x1, double h0, double h1, double g0, double g1){
  //computes x coord of point of intersection of tangents at x0 and x1
  //h0, h1 = heights at two points
  //g0, g1 = gradients at the two points
  if(fabs(g0 - g1) < EPS){
    cerr<<"Error in args passed to TangentIntersection: g0 = "<<g0<<", g1 = "<<g1<<endl;
    exit(1);
  }
  return (h1 - h0 - x1*g1 + x0*g0) / (g0 - g1);
}

double NewARS::UpperHull(double x, double x0, double h, double g){
  //computes value of upper hull at x
  //x0 = the point in the same interval as x
  //h, g = height and gradient at x0
  return h + (x-x0)*g;
}

double NewARS::ARS(const double* const args){

  //1: Initialise points, x. Typically the mode and a point either side
  double x[] = {-1.0, 0.0, 1.0};//dummy args

  //some code to find initial points goes here

  InitialisePoints(x, args);

  //?? no need to Initialise sampling function s
  //?? no need to initialise lower hull ??

  bool success = false;
  do{
  //Alternately
  //SAMPLE:
    success = Sample(args);

    if(!success){
  //UPDATE
      Update();
    }
  }
  while (!success);

  return NewPoint.abscissa;
}

void NewARS::test(){
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


int main() {
  NewARS Sampler;
  Sampler.Initialise(false, false, 0, 0, LogNormalDensity, LogNormalGradient);
  double args[2] = {0.0, 1.0};
  for(unsigned i = 0; i < 10; ++i)
  cout<<"Sampled point: "<<Sampler.ARS(args)<<endl<<endl;

  //Sampler.test(); 
}
