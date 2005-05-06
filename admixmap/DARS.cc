#include "DARS.h"

using namespace std;

// c *************************************************************
// c
// c    This program uses adaptative rejection algorithm for  
// c    log-concave distributions to generate from the  
// c	distribution of w in (-Infty,Infty)

// c 	Supply:  no:  Number of starting points
// c		lgth: Number of maxium points on the grid
// c	        alpha: vector of parameters
// c
// c	Externals: h  log-density
// c		   dh  derivative of the log-density
// c
// c		   x vector that contains the grid points
// c		   z vector of intersections of the tangent
// c		     lines
// c		   u=f(z)
// c
// c	Return
// c		   flag Counts the points used to get w	
// c		   w sample point


DARS::DARS()
{
   no = 3;
   loc = 0;
   lgth = 25;
   x0 = 0;
   f = new double[ lgth ];  
   df = new double[ lgth ];  
   x = new double[ lgth ];  
   u = new double[ lgth ];  
   z = new double[ lgth ];  
   psum = new double[ lgth ];
   MatrixArray_i null_MatrixArray_i(1);
   MatrixArray_d null_MatrixArray_d(1);
   data_i = null_MatrixArray_i;
   data_d = null_MatrixArray_d;
}

DARS::DARS( int inLeftFlag, int inRightFlag, double innewnum,
            //const Vector_d &inparameters,
	    const double inparameters[],int size,
            double (*funct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
            double (*dfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
            double (*ddfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
            const MatrixArray_i &integer_data, const MatrixArray_d &double_data )
{
   no = 3;
   loc = 0;
   lgth = 25;
   x0 = 0;
   //parameters = inparameters;
   parameters.SetNumberOfElements(size);
   for(int i=0;i<size;++i)parameters(i)=inparameters[i];
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

void DARS::
SetParameters( int inLeftFlag, int inRightFlag, double innewnum,
               //const Vector_d &inparameters,
	       const double inparameters[],int size,
               double (*funct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
               double (*dfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
               double (*ddfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
               const MatrixArray_i &integer_data, const MatrixArray_d &double_data )
{
  //parameters = inparameters;
  parameters.SetNumberOfElements(size);
  for(int i=0;i<size;++i)parameters(i)=inparameters[i];
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

//void DARS::UpdateParameters( const Vector_d &inparameters )
void DARS::UpdateParameters( const double inparameters[], int size )
{
  //parameters = inparameters;
  for(int i=0;i<size;++i)parameters(i)=inparameters[i];
}

void DARS::UpdateIntegerData( const MatrixArray_i &indata )
{
   data_i = indata;
}

void DARS::UpdateDoubleData( const MatrixArray_d &indata )
{
   data_d = indata;
}

void DARS::BeginModeSearch( double innewnum )
{
   newnum = innewnum;
}

double DARS::Sample()
{
   double w, dfa = 1, dfb = -1;

   if( !LeftFlag )
      dfa = (*dfunction)( parameters, data_i, data_d, x0 + 0.00000001 );
   if( !RightFlag )
      dfb = (*dfunction)( parameters, data_i, data_d, x1 - 0.00000001 );

// Check gradient is positive at minimum value and negative at maximum value
   if( dfa > 0 && dfb < 0 ){
      
// Mode search
      if( !LeftFlag && !RightFlag ){
         SimpleModeSearch( x0, x1 );
      }
      else{
         NewtonRaphson();
      }
      while( isinf( (*function)( parameters, data_i, data_d, x[0] ) ) )
         x[0] = ( x[0] + x[1] ) / 2;
      while( isinf( (*function)( parameters, data_i, data_d, x[2] ) ) )
         x[2] = ( x[2] + x[1] ) / 2;
   }
   else if( dfb > 0 ){
      newnum = x1 - 0.00000001;
      x[2] = newnum;
      if( !LeftFlag ){
         x[0] = x0 + 0.00000001;
         x[1] = (x[0] + x[2]) / 2;
      }
      else{
         x[1] = x[2] - 2/dfb;
         x[0] = x[1] - 2/dfb;
      }
   }
   else{
      newnum = x0 + 0.00000001;
      x[0] = newnum;
      if( !RightFlag ){
         x[2] = x1 - 0.00000001;
         x[1] = (x[0] + x[2]) / 2;
      }
      else{
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
   double aux1 = 0.0;
   double bux,max,un,w,temp,newdf;
         
   n = no;
   for( int i = 0; i < lgth; i++ )
   {
      psum[i] = 0;
      u[i] = 0;
      z[i] = 0;
   }

   max=(*function)(parameters, data_i, data_d, newnum );

   for( int i = 0; i < n; i++ )
   {
      aux=x[i];
      f[i]=(*function)(parameters, data_i, data_d, aux)-max;
      df[i]=(*dfunction)(parameters, data_i, data_d, aux);
   }


   loc = 0;

// ** Routine starts **

   do
   {
      consz();
      conspsum();
      
      un=myrand();
      
      w=un*psum[n-1];
      
      i = 0;
      
      while( i < n && w > psum[i] )
         i++;
      
      loc = i;
         
//  ** Generate w and height of w on the envelope distribution  **

      if( loc == 0 )
      {
         if( LeftFlag == 0 ) // Ok for sampling from x0
            aux = w*df[0]+fannyexp(df[0]*(x0-x[0])+f[0]);
         else if( LeftFlag == 1 ) // Sample from -Inf
            aux=w*df[0];
         w = x[0]+(log(aux)-f[0])/df[0];
         if( isnan(w) ){
            cout << "Nan at loc = 0\n";
         }
      }
      else if( fabs(df[loc]) < 0.0001 )
      {
         w -= psum[loc-1];
         aux = fannyexp(f[loc]);
         w = z[loc-1]+w/aux;
         if( isnan(w) ){
            cout << "Nan at df[loc] < 0.0001 = 0\n";
         }
      }
      else   
      {
         w -= psum[loc-1];
         aux = w*df[loc]+fannyexp(u[loc-1]);         
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
         aux1 = (*function)(parameters, data_i, data_d, w)-max;
         bux = fannyexp(aux1);
         if ( log(un) <= aux1 - log(aux) )
            SampleFlag = 1;
      }

// ** Prepare new envelope **

      newdf = (*dfunction)(parameters, data_i, data_d, w);
      if( SampleFlag == 0 ){
         if( !isinf(aux1) && !isinf(newdf) && !isnan(aux1) && !isnan(newdf) ){
         n++;
         
         if( n >= lgth )
         {
            cout << "warning" << n << endl;
            cout << x[0] << " " << x[n-1] << endl;
            exit(1);
         }
         
         for( i = n - 1; i >= loc + 1; i-- )
         {
            x[i] = x[i-1];
            f[i] = f[i-1];
            df[i] = df[i-1];
         }
         x[loc] = w;
         f[loc] = aux1;
         df[loc] = newdf;
         loc--;
         }
         else{
            cout << "bollox\n";
         }
      }
   }while( SampleFlag == 0 );
   
   flag += n-no+1;

   return ( w );
}

void DARS::consz()
{
   int i,iloc,iloc1;
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
         u[i] = df[i]*(z[i]-x[i])+f[i];
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
      z[n-1] = x1;
      u[n-1] = df[n-1]*(z[n-1]-x[n-1])+f[n-1];
   }
}

void DARS::conspsum()
{
   int i,iloc,iloc1;
   double aux1,aux2,saux1;

   double aux = 0.0;
   double saux = 0.0;
 
   aux1 = psum[0];
   if( loc == -1 )
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
   else if( loc == 0 )
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

   for( i = iloc; i  < n - RightFlag; i++ )
   {
      if( (n == no) || (i <= iloc1) )
      {
         aux = aux1;
         aux1 = psum[i];
         if ( fabs(df[i]) < 0.0001 )
         {
            saux1 = (z[i]-z[i-1])*fannyexp(f[i]);
            saux += saux1;
         }
         else
         {
            saux1 = (fannyexp(u[i])-fannyexp(u[i-1]))/df[i];
            saux += saux1;
         }
         psum[i] = saux;
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
      else
      {
         aux2 = psum[i];
         saux = aux1-aux+psum[i-1];
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
   }while( ( (dfold * dfnew) > 0 ) && fabs(dfnew) > 0.01 );

   if( fabs(dfnew) > 0.01 ){
      if( oldnum < newnum ){
         a = oldnum;
         b = newnum;
      }
      else{
         a = newnum;
         b = oldnum;
      }
      
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
