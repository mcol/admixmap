#include "ModeFinder.h"
#include <math.h>
#include <iostream>

double ModeFinder::FindModeByNewtonRaphson(double init, double* ddf, const void* const args, 
					   double (*function)(double, const void* const), 
					   double (*dfunction)(double, const void* const),
					   double (*ddfunction)(double, const void* const))
{
  //init is an initial guess, args are arguments to loglikelihood and derivative functions
  //returns mode as estimate of mode of f
  double step = 0.0, f = 0.0, flast = 0.0;
  double df = (*dfunction)( init, args );//1st derivative at initial value
  *ddf = (*ddfunction)( init, args );//2nd derivative " " 

  double newnum = init;
  while( fabs(df) > 0.001 ){
    flast = (*function)(newnum, args);//loglikelihood at current value
    step = -df / *ddf;

    f = (*function)(step+newnum, args);//loglikelihood at new value
    //std::cout << "df = " << df << " ddf= " << *ddf <<  " f = " << f << " " << f-flast << " step = " << step <<  " newnum = " << newnum << std::endl;
     if(f >= flast)//if NR step increases loglikelihood
     newnum += step;//use it
     else{//step is too large and we've overshot the mode so try a smaller step
       do{
 	step *= 0.5;
// 	if(fabs(step) < 0.000001);
// 	step *= 100.0;
 	f = (*function)(newnum+step, args);
 	}
 	while(f < flast);//repeat until loglikelihood increases
       newnum += step;
     }

    df = (*dfunction)( newnum, args );//update 1st derivative
    *ddf = (*ddfunction)( newnum, args );//update 2nd derivative
  }
  //std::cout << "mode found at " << newnum << std::endl;
  return newnum;
}

double ModeFinder::SimpleModeSearch( const void* const args,
				     double (*gradient)(double, const void* const), double (*secondDeriv)(double, const void* const),
				     bool hasLowerBound, bool hasUpperBound, const double LowerBound, const double UpperBound ){
  double mode = 0.0;
 if( !hasLowerBound || !hasUpperBound )
   mode = UnboundedModeSearch(args, gradient, secondDeriv, hasLowerBound, hasUpperBound, LowerBound, UpperBound );

 else BoundedModeSearch(LowerBound, UpperBound, args, gradient);
  return mode;
}

///searches for mode between a and b, without evaluating gradient at a or b 
double ModeFinder::BoundedModeSearch( const double a, const double b, const void* const args,
				     double (*gradient)(double, const void* const) )

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
  // step 2: iterate using average 2nd derivative for approximate Newton-Raphson update 
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

/**
   Uses modified Newton-Raphson (step size * 2) to find two points
   either side of the mode. Given these two points use most simple and
   robust mode search SimpleModeSearch(). Two points chosen such that
   they are 3 * second derivative from the mode.
*/
double ModeFinder::UnboundedModeSearch(const void* const args, double (*gradient)(double, const void* const),
				       double (*secondDeriv)(double, const void* const), bool hasLowerBound, bool hasUpperBound, 
				       const double LowerBound, const double UpperBound )
{
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
    
    mode = BoundedModeSearch( a, b, args, gradient);

  }//end if(fabs(dfnew > 0.01)
  else{
    mode = newnum;
  }
  return mode;
}
