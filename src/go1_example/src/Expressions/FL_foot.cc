/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:28 GMT-05:00
 */

#include "math2mat.hpp"
#include "mdefs.hpp"

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t155;
  double t249;
  double t341;
  double t403;
  double t475;
  double t185;
  double t648;
  double t115;
  double t711;
  double t719;
  double t907;
  double t952;
  double t1016;
  double t1045;
  double t1161;
  double t1223;
  double t1229;
  double t1269;
  double t1357;
  double t1455;
  double t1463;
  double t130;
  double t133;
  double t1677;
  double t414;
  double t448;
  double t1759;
  double t1801;
  double t1807;
  double t1808;
  double t1148;
  double t1174;
  double t1206;
  double t1726;
  double t1772;
  double t1779;
  double t1922;
  double t1924;
  double t1932;
  double t1309;
  double t1328;
  double t1352;
  double t2024;
  double t2026;
  double t2029;
  double t2063;
  double t2077;
  double t2080;
  double t2398;
  double t2410;
  double t2414;
  double t2304;
  double t2310;
  double t2338;
  double t2496;
  double t2513;
  double t2517;
  double t2545;
  double t2578;
  double t2581;
  double t2604;
  double t2606;
  double t2617;
  t155 = Cos(var1[4]);
  t249 = Sin(var1[9]);
  t341 = Sin(var1[4]);
  t403 = Cos(var1[9]);
  t475 = Sin(var1[5]);
  t185 = Cos(var1[5]);
  t648 = Sin(var1[10]);
  t115 = Cos(var1[10]);
  t711 = t403*t341;
  t719 = t155*t249*t475;
  t907 = t711 + t719;
  t952 = Cos(var1[11]);
  t1016 = -1.*t952;
  t1045 = 1. + t1016;
  t1161 = Sin(var1[11]);
  t1223 = t155*t185*t648;
  t1229 = t115*t907;
  t1269 = t1223 + t1229;
  t1357 = t115*t155*t185;
  t1455 = -1.*t648*t907;
  t1463 = t1357 + t1455;
  t130 = -1.*t115;
  t133 = 1. + t130;
  t1677 = Sin(var1[3]);
  t414 = -1.*t403;
  t448 = 1. + t414;
  t1759 = Cos(var1[3]);
  t1801 = t1759*t185;
  t1807 = -1.*t1677*t341*t475;
  t1808 = t1801 + t1807;
  t1148 = -0.35*t1045;
  t1174 = 0.3455*t1161;
  t1206 = t1148 + t1174;
  t1726 = t185*t1677*t341;
  t1772 = t1759*t475;
  t1779 = t1726 + t1772;
  t1922 = -1.*t403*t155*t1677;
  t1924 = -1.*t249*t1808;
  t1932 = t1922 + t1924;
  t1309 = 0.3455*t1045;
  t1328 = 0.35*t1161;
  t1352 = t1309 + t1328;
  t2024 = t648*t1779;
  t2026 = t115*t1932;
  t2029 = t2024 + t2026;
  t2063 = t115*t1779;
  t2077 = -1.*t648*t1932;
  t2080 = t2063 + t2077;
  t2398 = t185*t1677;
  t2410 = t1759*t341*t475;
  t2414 = t2398 + t2410;
  t2304 = -1.*t1759*t185*t341;
  t2310 = t1677*t475;
  t2338 = t2304 + t2310;
  t2496 = t403*t1759*t155;
  t2513 = -1.*t249*t2414;
  t2517 = t2496 + t2513;
  t2545 = t648*t2338;
  t2578 = t115*t2517;
  t2581 = t2545 + t2578;
  t2604 = t115*t2338;
  t2606 = -1.*t648*t2517;
  t2617 = t2604 + t2606;
  p_output1[0]=t1206*t1269 + t1352*t1463 + 0.3455*t133*t155*t185 - 0.072*t249*t341 - 0.072*t155*t448*t475 + 0.19875*(t249*t341 - 1.*t155*t403*t475) + 0.3455*t648*t907 - 0.7*(t1161*t1463 + t1269*t952) + 0.3455*(-1.*t1161*t1269 + t1463*t952) + var1[0];
  p_output1[1]=0.3455*t133*t1779 + t1206*t2029 + t1352*t2080 + 0.072*t155*t1677*t249 + 0.19875*(-1.*t155*t1677*t249 + t1808*t403) + 0.072*t1808*t448 + 0.3455*t1932*t648 - 0.7*(t1161*t2080 + t2029*t952) + 0.3455*(-1.*t1161*t2029 + t2080*t952) + var1[1];
  p_output1[2]=0.3455*t133*t2338 - 0.072*t155*t1759*t249 + t1206*t2581 + t1352*t2617 + 0.19875*(t155*t1759*t249 + t2414*t403) + 0.072*t2414*t448 + 0.3455*t2517*t648 - 0.7*(t1161*t2617 + t2581*t952) + 0.3455*(-1.*t1161*t2581 + t2617*t952) + var1[2];
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 18 && ncols == 1) && 
      !(mrows == 1 && ncols == 18))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 3, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "FL_foot.hh"

namespace SymFunction
{

void FL_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
