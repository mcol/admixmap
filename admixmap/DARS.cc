/** 
 *   ADMIXMAP
 *   DARS.cc
 *   This class implements an adaptative rejection algorithm for  
 *   log-concave distributions to generate from the  
 *   distribution of w in (-Infty,Infty)
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
#include "DARS.h"

using namespace std;


/*
  Supply:  no:  Number of starting points
  lgth: Number of maxium points on the grid
  alpha: vector of parameters
  
  Externals: h  log-density
  dh  derivative of the log-density
  
  x vector of points on x axis at which gradient is evaluated
  z vector of points on x axis at which i th tangent intersects i+1 th tangent
  lines
  u  vector of heights at which i th tangent intersects i+1th tangent   

  LeftFlag and RightFlag are coded as 0 if there is a left truncation point, 1 otherwise
  left and right truncation points are passed separately as inx0 and inx1 

  Returns
  flag Counts the points used to get w	
  w sample point
*/

// should have the following functions:-
// SetStartingValues: sets 3 starting values and initializes upper and lower hulls
// SampleFromUpperHull: samples from distribution defined by upper hull
// HeightUpperHull: evaluates height of upper hull at position x
// HeightLowerHull: evaluates height of lower hull at position x
// UpdateHulls: updates upper and lower hulls with new tangent and chord at position x
//    this should include rearranging the x values in ascending order 
//    and recalculating areas under curves formed by tangents   
// algorithm samples x from SampleUpperHull and w from unif(0,1) 
// squeezing test calls functions HeightUpperHull and HeightLowerHull
// rejection test evaluates log density at x
// if rejected by both tests, call UpdateHulls and repeat  
// how is distribution defined by upper hull normalized?  
// should use STL vectors to allow x values, heights and gradients to be inserted in sequence

DARS::DARS()
{
   no = 3;
   loc = 0;
   lgth = 25;
   x0 = 0;
   f = new double[ lgth ];  
   df = new double[ lgth ];  // vector of gradients at x[]
   x = new double[ lgth ];  // vector of parameter values at which gradient is evaluated
   u = new double[ lgth ];  // u[i] is log density at z[i] 
   z = new double[ lgth ];  // vector of parameter values at which tangents intersect
   psum = new double[ lgth ];
   data_i = 0;
   data_d = 0;
   parameters = 0;
}

DARS::DARS( int inLeftFlag, int inRightFlag, double innewnum,
	    const double inparameters[],
            double (*funct)(const double*, const int *, const double*, double),
            double (*dfunct)(const double*, const int *, const double*, double),
            double (*ddfunct)(const double*, const int *, const double*, double),
            const int *integer_data, const double *double_data )
{
  no = 3; // num starting points
  loc = 0;
  lgth = 25; // max num points 
  x0 = 0;
  parameters = inparameters;
  data_i =  integer_data;
  data_d =  double_data;
  function = funct;
  dfunction = dfunct;
  ddfunction = ddfunct;
  f = new double[ lgth ];
  df = new double[ lgth ];  
  x = new double[ lgth ];  
  u = new double[ lgth ];  
  z = new double[ lgth ];
  psum = new double[ lgth ];
  LeftFlag = inLeftFlag;
  RightFlag = inRightFlag;
  newnum = innewnum;
}

DARS::~DARS()
{
   delete [] psum;
   delete [] x;
   delete [] f;
   delete [] df;
   delete [] z;
   delete [] u;
}

void DARS::SetParameters( int inLeftFlag, int inRightFlag, double innewnum,
	       const double inparameters[],
               double (*funct)(const double*, const int *, const double*, double),
               double (*dfunct)(const double*, const int *, const double*, double),
               double (*ddfunct)(const double*, const int *, const double*, double),
               const int *integer_data, const double *double_data )
{
  parameters = inparameters;
  data_i = integer_data;
  data_d = double_data;
  function = funct;
  dfunction = dfunct;
  ddfunction = ddfunct;
  LeftFlag = inLeftFlag;
  RightFlag = inRightFlag;
  newnum = innewnum;
}

void DARS::SetLeftTruncation( double inx0 )
{
   x0 = inx0;
}

void DARS::SetRightTruncation( double inx0 )
{
   x1 = inx0;
}

void DARS::UpdateParameters( const double inparameters[])
{
  parameters = inparameters;
}

void DARS::UpdateIntegerData( const int *indata )
{
   data_i = indata;
}

void DARS::UpdateDoubleData( const double *indata )
{
   data_d = indata;
}

void DARS::BeginModeSearch( double innewnum )
{
   newnum = innewnum;
}

double DARS::Sample()
  // function should have argument allowing user to pass mode if known 
{
  double w, dfa = 1, dfb = -1;
  
  if( !LeftFlag ) // if bounded below, evaluate gradient just above minimum value 
    dfa = (*dfunction)( parameters, data_i, data_d, x0 + 0.00000001 );
  if( !RightFlag ) // if bounded above, evaluate gradient just below max value
    dfb = (*dfunction)( parameters, data_i, data_d, x1 - 0.00000001 );
  
  // Check gradient is positive at minimum value and negative at maximum value
  if( dfa > 0 && dfb < 0 ){
    
    // Mode search
    if( !LeftFlag && !RightFlag ){ // if bounded below and above 
      SimpleModeSearch( x0, x1 );
    }
    else{ // if unbounded below or above
      NewtonRaphson(); // assigns x[1] as mode, x[0], x[2] as +/- 3 times 2nd deriv
    }
    
    // if log density infinite at x[0], assign x[0] as mean of x[0] and x[1]
    while( isinf( (*function)( parameters, data_i, data_d, x[0] ) ) )
      x[0] = ( x[0] + x[1] ) / 2;
    // x[2] similarlyx
    while( isinf( (*function)( parameters, data_i, data_d, x[2] ) ) )
      x[2] = ( x[2] + x[1] ) / 2;
  }
  else if( dfb > 0 ){ // if gradient positive at maximum value
    newnum = x1 - 0.00000001;
    x[2] = newnum;
    if( !LeftFlag ){ // if bounded below
      x[0] = x0 + 0.00000001;
      x[1] = (x[0] + x[2]) / 2; // assign x[1] as mean of x[0] and x[2]
    }
    else{ // use gradient at max to assign x[0] and x[2]
      x[1] = x[2] - 2/dfb;
      x[0] = x[1] - 2/dfb;
    }
  }
  else{ // if gradient negative at minimum value, 
    newnum = x0 + 0.00000001;
    x[0] = newnum;
    if( !RightFlag ){ // assign x[1] as mean of x[0] and x[2]
      x[2] = x1 - 0.00000001;
      x[1] = (x[0] + x[2]) / 2;
    }
    else{ // use gradient at min
      x[1] = x[0] - 2/dfa; 
      x[2] = x[1] - 2/dfa;
    }
  }
  
  w = SampleUsingARS();
  return( w );
}

double DARS::SampleUsingARS()
{
  int flag = 0, SampleFlag = 0, i;
  double aux = 0.0;
  double aux1 = 0.0; // log density at x, minus max 
  double bux,max,un,w,temp,newdf; // bux is exp(aux1)
  // max is log density at newnum

  n = no; // number of points at which log density and gradient are evaluated, initially set to num starting points
  for( int i = 0; i < lgth; i++ ) // loop from 0 to max num grid points
    {
      psum[i] = 0; // i th element is sum of areas under tangents 0 to i. 
      u[i] = 0;
      z[i] = 0; // 
    }
  
  max=(*function)(parameters, data_i, data_d, newnum ); // log density at newnum
  // loop over n grid points to calculate log density f[i] and gradient df[i] 
  for( int i = 0; i < n; i++ )
    {
      aux=x[i]; // here aux is assigned value of x[i]
      //  elsewhere it is used for height of log density at x[i]
      f[i]=(*function)(parameters, data_i, data_d, aux)-max;
      df[i]=(*dfunction)(parameters, data_i, data_d, aux);
    }
  
  loc = 0;
  // ** Routine starts **
  do
    {
      consz();  // update vectors specifying tangents: horizontal positions z and heights u at intersections 
      conspsum();  // update vector of cumulative sums of areas under density curves formed by tangents to log density
      
      // draw w as uniform between 0 and psum[n-1]: here w is a value from the cumulative distribution function 
      un=myrand();
      w=un*psum[n-1]; 
      // assign loc as integer between 0 and n that specifies which tangent w is under
      i = 0;  
      while( i < n && w > psum[i] )
	i++;
      loc = i;
      
      //  ** convert w to position on x axis, and calculate aux: height of envelope at position w   **
      if( loc == 0 ) // if no tangents intersect to left of w
	{
	  if( LeftFlag == 0 ) // if lower bound
            aux = w*df[0] + fannyexp(df[0]*(x0 - x[0]) + f[0]); // height of envelope at lower bound
	  else if( LeftFlag == 1 ) // if no lower bound 
            aux=w*df[0]; // height of envelope at w=0    

	  w = x[0]+(log(aux)-f[0])/df[0]; // w recoded as parameter value
	  if( isnan(w) ){
            cout << "Nan at loc = 0\n";
	  }
	}
      else if( fabs(df[loc]) < 0.0001 ) // if gradient flat at x[loc]
	{
	  w -= psum[loc-1]; 
	  aux = fannyexp(f[loc]);  // height of density at x[loc]  
	  w = z[loc-1]+w/aux;  // w recoded as parameter value 
	  if( isnan(w) ){
            cout << "Nan at df[loc] < 0.0001 = 0\n";
	  }
	}
      else  // at least one tangent to left of w, and gradient not flat at x[loc]
	{
	  w -= psum[loc-1];
	  aux = w*df[loc]+fannyexp(u[loc-1]);  // density at z[loc-1 + height of tangent at w=0       
	  w = x[loc]+(log(aux)-f[loc])/df[loc];
	  if( isnan(w) ){
            cout << "Nan at else = 0\n";
	  }
	}

      //** Prepare squeezing pre-test **
      if ( w >= x[loc] && loc == n - 1 )
	{
	  bux=0.0;
	  loc++;
	}
      else if( w >= x[loc] )
	{
	  temp=(f[loc+1]*(w-x[loc])) / (x[loc+1]-x[loc]);
	  bux=fannyexp(temp);
	  bux *= fannyexp(f[loc]*(x[loc+1]-w));
	  loc++;
	}
      else if( loc == 0 )
	{
	  bux=0.0;
	}
      else
	bux = fannyexp((f[loc-1]*(x[loc]-w)+f[loc] * (w-x[loc-1]))/(x[loc]-x[loc-1]));
      
      // ** Rejection step by squeezing  **
      un = myrand();
      if( un <= (bux/aux) )
	SampleFlag = 1;
      else
	{
	  // ** Rejection step  **
	  aux1 = (*function)(parameters, data_i, data_d, w) - max; // log density at w, minus max
	  bux = fannyexp(aux1);  // density at w
	  if ( log(un) <= aux1 - log(aux) )  // aux is height of envelope at w
            SampleFlag = 1; // exit loop
	}
      
      // ** Prepare new envelope **
      newdf = (*dfunction)(parameters, data_i, data_d, w); // log density at parameter value w 
      if( SampleFlag == 0 ){
	if( !isinf(aux1) && !isinf(newdf) && !isnan(aux1) && !isnan(newdf) ){ // if aux1 and newdf not inf or nan 
	  n++;
	  
	  if( n >= lgth ) // failure to accept after max number of grid points
	    {
	      cout << "Adaptive rejection sampler: failure to accept values " << n << endl;
	      cout << x[0] << " " << x[n-1] << endl;
	      exit(1);
	    }
	  
	  for( i = n - 1; i >= loc + 1; i-- )
	    { // loop from i=n-1 to i=loc+1
	      x[i] = x[i-1]; 
	      f[i] = f[i-1];
	      df[i] = df[i-1];
	    }
	  x[loc] = w; // update vector of grid points with w
	  f[loc] = aux1; // update vector of log densities at x
	  df[loc] = newdf; // update vector of gradients at x
	  loc--; // decrement loc - this could set loc to -1 if no tangents intersect to left of w
	}
	else{ // aux1 or newdf inf or nan
	  cout << "Adaptive rejection sampler: log density or gradient infinite or NaN\n";
	}
      }
   }while( SampleFlag == 0 );
   
   flag += n-no+1;
   
   return ( w );
}

void DARS::consz() // constructs envelope of tangents 
  // updates scalars aux, bux, aux1, aux2, arrays z, u
{
  int i,iloc,iloc1;
  // same names as global vars
  double aux,bux;
  double aux1 = 0.0;
  double bux1 = 0.0;
  
  if( loc == -1 )
    {
      iloc=0;
      iloc1=0;
    }
  else
    {
      iloc=loc;
      iloc1=loc+1;
    }
  
  for( i = iloc; i < n - 1; i++ )
    {
      if( (n == no) || (i <= iloc1) )
	{
	  aux1 = z[i];
	  bux1 = u[i];
	  z[i] = (f[i+1]-f[i]-x[i+1]*df[i+1]+x[i]*df[i])/ (df[i]-df[i+1]); 
	  //horizontal position of intersection of i th tangent with (i+1) th tangent  
	  u[i] = df[i]*(z[i]-x[i])+f[i]; // height of i th tangent at horizontal position z[i]
	}
      else
	{
	  aux = z[i];
	  bux = u[i];
	  z[i] = aux1;
	  u[i] = bux1;
	  aux1 = aux;
	  bux1 = bux;
	}
    }
  if( !RightFlag ){
    z[n-1] = x1; // horizontal position of last intersection is set to upper bound
    u[n-1] = df[n-1]*(z[n-1]-x[n-1])+f[n-1];
  }
}

void DARS::conspsum()
{
  // updates psum: array of cumulative sums of areas under curves formed by tangents to log density 
  // should take as const arguments: x[], z[], u[], df[], x0. x1, loc
  // new x value is under the loc th tangent
  // in case of error, should output all arguments 
  // 
  int i,iloc,iloc1;
  double aux1,aux2,saux1;
  
  double aux = 0.0; // stores last value of aux1
  double saux = 0.0; // when calculated, this is assigned to psum[i]
  
  aux1 = psum[0]; // stores last value of aux2

  if( loc == -1 ) // no tangents intersect to left of w, and decrement after rejection step
    {
      if( LeftFlag == 0 ) // Ok for sampling from x0
	saux = fannyexp(u[0])*(1 - fannyexp((x0-z[0])*df[0]))/df[0];
      else if( LeftFlag == 1 ) // Sample from -Inf
	saux = fannyexp(u[0])/df[0];
      psum[0] = saux;
      iloc = 1;
      iloc1 = 1;
      if( isinf(saux) ){
	cout << "suax = inf at loc == -1" << endl;
	exit(0);
      }
    }
  else if( loc == 0 ) // no tangents intersect to left of w, or loc=1 before decrement
    {
      if( LeftFlag == 0 ) // Ok for sampling from x0
	saux = fannyexp(u[0])*(1. - fannyexp((x0-z[0])*df[0]))/df[0];
      else if( LeftFlag == 1 ) // Sample from -Inf
	saux = fannyexp(u[0])/df[0];
      psum[0] = saux;
      iloc = 1;
      iloc1 = 2;
      if( isinf(saux) ){
	cout << "suax = inf at loc == 0" << endl;
	exit(0);
      }
    }
  else  
    {
      saux = psum[loc-1];
      iloc = loc;
      iloc1 = loc+2;
      if( isinf(saux) ){
	cout << "suax = inf at loc != -1, 0." << endl;
	exit(0);
      }
    }
  
  for( i = iloc; i  < n - RightFlag; i++ ) // loops from iloc to n-1 if no upper bound, to n-2 if no upper bound 
    {
      if( (n == no) || (i <= iloc1) ) 
	{
	  aux = aux1;
	  aux1 = psum[i];
	  if ( fabs(df[i]) < 0.0001 ) // if gradient is flat at i th x value
	    {
	      saux1 = (z[i]-z[i-1])*fannyexp(f[i]);
	      saux += saux1;
	    }
	  else
	    {
	      saux1 = (fannyexp(u[i])-fannyexp(u[i-1]))/df[i]; // positive if numerator has same sign as denominator 
	      saux += saux1;
	    }
	  psum[i] = saux; // assign i th element of psum 
	  if( psum[i] < 0 ){ 
            cout << "psum < 0 : 1" << endl;
            cout << n << endl;
            for( int ii = 0; ii < n; ii++ )
	      cout << x[ii] << " ";
            cout << endl;
            for( int ii = 0; ii < n; ii++ )
	      cout << f[ii] << " ";
            cout << endl;
            for( int ii = 0; ii < n; ii++ )
	      cout << df[ii] << " ";
            cout << endl;
	  }
	}
      else // n not equal to num starting points and i > loc + 2
	{
	  aux2 = psum[i];
	  saux = aux1 - aux + psum[i-1];
	  psum[i] = saux;
	  if( psum[i] < 0 ){
            cout << "psum < 0 : 2" << endl;
            cout << n << endl;
            for( int ii = 0; ii < n; ii++ )
	      cout << x[ii] << " ";
            cout << endl;
            for( int ii = 0; ii < n; ii++ )
	      cout << f[ii] << " ";
            cout << endl;
            for( int ii = 0; ii < n; ii++ )
	      cout << df[ii] << " ";
            cout << endl;
	  }
	  else if( isinf(psum[i]) || isinf(psum[i]) ){
            cout << "Error: psum[i] = inf" << endl;
            exit(0);
	  }
	  aux = aux1;
	  aux1 = aux2;
	}
    }
  // Change in case you are not sampling from Infty
  if( RightFlag )
    psum[n-1] = saux - fannyexp(u[n-2])/df[n-1];
  if( psum[n-1] < 0 ){
    cout << "psum < 0 : 3" << endl;
    cout << n << endl;
    for( int ii = 0; ii < n; ii++ )
      cout << x[ii] << " ";
    cout << endl;
    for( int ii = 0; ii < n; ii++ )
      cout << f[ii] << " ";
    cout << endl;
    for( int ii = 0; ii < n; ii++ )
      cout << df[ii] << " ";
    cout << endl;
  }
  else if( isinf(psum[n-1]) || isinf(psum[n-1]) ){
    cout << "Error: psum[n-1] = inf" << endl;
    exit(0);
  }
}

void DARS::SimpleModeSearch( double aa, double bb )
{
   double a, b, num, df2, df1, dfa, dfb, ddf;
   int count = 0;
   a = aa + 0.00000001;
   b = bb - 0.00000001;

   newnum = ( a + b ) / 2;
   do{
      count++;
      df1 = (*dfunction)(parameters, data_i, data_d, newnum);
      if( df1 > 0.0 ){
         num = (newnum + b) / 2;
         df2 = (*dfunction)(parameters, data_i, data_d, num);
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
         df2 = (*dfunction)(parameters, data_i, data_d, num);
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
   dfa = (*dfunction)(parameters, data_i, data_d, a);
   dfb = (*dfunction)(parameters, data_i, data_d, b);
   x[1] = newnum;

   if( dfa > 1.0 )
      x[0] = a;
   else{
      ddf = (*ddfunction)( parameters, data_i, data_d, newnum );
      x[0] = x[1] + 3.0 / ddf;
      if( LeftFlag == 0 && x[0] < x0 )
         x[0] = x0 + 0.00000001;
   }

   if( dfb < -1.0 )
      x[2] = b;
   else{
      ddf = (*ddfunction)( parameters, data_i, data_d, newnum );
      x[2] = x[1] - 3.0 / ddf;
      if( RightFlag == 0 && x[2] > x1 )
         x[2] = x1 - 0.00000001;
   }
}

void DARS::NewtonRaphson()
{
// Using modified Newton-Raphson (step size * 2) to find two points
// either side of the mode. Given these two points use most simple and
// robust mode search SimpleModeSearch(). Two points chosen such that
// they are 3 * second derivative from the mode.
   double oldnum, step, dfnew, dfold, ddf, a, b;

   do{
      oldnum = newnum;
      ddf = (*ddfunction)( parameters, data_i, data_d, oldnum );
      dfold = (*dfunction)( parameters, data_i, data_d, oldnum );
      step = -dfold / ddf;
      newnum += 2 * step;
      if( !LeftFlag && newnum < x0 )
         newnum = x0 + 0.00000001;
      if( !RightFlag && newnum > x1 )
         newnum = x1 - 0.00000001;
      dfnew = (*dfunction)( parameters, data_i, data_d, newnum );
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
      // if a below lower bound or b below upper bound, reset to just within bounds 
      if( LeftFlag == 0 && a < x0 )
         a = x0 + 0.00000001;
      else if( RightFlag == 0 && b > x1 )
         b = x1 - 0.00000001;
      
      SimpleModeSearch( a, b );
   }
   else{
      x[1] = newnum;
      x[0] = x[1] + 3.0 / ddf;
      x[2] = x[1] - 3.0 / ddf;
      if( LeftFlag == 0 && x[0] < x0 )
         x[0] = x0 + 0.00000001;
      else if( RightFlag == 0 && x[2] > x1 )
         x[2] = x1 - 0.00000001;
   }
}


double
DARS::fannyexp( double x )
{
   double y;
   if( x > -700 )
      y = exp(x);
   else
      y = 0;
   return( y );
}
