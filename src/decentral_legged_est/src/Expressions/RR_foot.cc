/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:30 GMT-05:00
 */

#include "math2mat.hpp"
#include "mdefs.hpp"

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t1503;
  double t1660;
  double t1707;
  double t1779;
  double t1904;
  double t1514;
  double t2094;
  double t1049;
  double t2179;
  double t2200;
  double t2221;
  double t2338;
  double t2364;
  double t2398;
  double t2422;
  double t2537;
  double t2578;
  double t2600;
  double t2727;
  double t2730;
  double t2748;
  double t1229;
  double t1466;
  double t2956;
  double t1832;
  double t1842;
  double t2996;
  double t3110;
  double t3114;
  double t3122;
  double t2410;
  double t2478;
  double t2496;
  double t2992;
  double t3036;
  double t3038;
  double t3255;
  double t3260;
  double t3262;
  double t2662;
  double t2705;
  double t2719;
  double t3268;
  double t3280;
  double t3289;
  double t3322;
  double t3324;
  double t3330;
  double t3482;
  double t3495;
  double t3497;
  double t3453;
  double t3457;
  double t3458;
  double t3520;
  double t3521;
  double t3527;
  double t3539;
  double t3657;
  double t3660;
  double t3726;
  double t3747;
  double t3749;
  t1503 = Cos(var1[4]);
  t1660 = Sin(var1[12]);
  t1707 = Sin(var1[4]);
  t1779 = Cos(var1[12]);
  t1904 = Sin(var1[5]);
  t1514 = Cos(var1[5]);
  t2094 = Sin(var1[13]);
  t1049 = Cos(var1[13]);
  t2179 = t1779*t1707;
  t2200 = t1503*t1660*t1904;
  t2221 = t2179 + t2200;
  t2338 = Cos(var1[14]);
  t2364 = -1.*t2338;
  t2398 = 1. + t2364;
  t2422 = Sin(var1[14]);
  t2537 = t1503*t1514*t2094;
  t2578 = t1049*t2221;
  t2600 = t2537 + t2578;
  t2727 = t1049*t1503*t1514;
  t2730 = -1.*t2094*t2221;
  t2748 = t2727 + t2730;
  t1229 = -1.*t1049;
  t1466 = 1. + t1229;
  t2956 = Sin(var1[3]);
  t1832 = -1.*t1779;
  t1842 = 1. + t1832;
  t2996 = Cos(var1[3]);
  t3110 = t2996*t1514;
  t3114 = -1.*t2956*t1707*t1904;
  t3122 = t3110 + t3114;
  t2410 = -0.35*t2398;
  t2478 = -0.3455*t2422;
  t2496 = t2410 + t2478;
  t2992 = t1514*t2956*t1707;
  t3036 = t2996*t1904;
  t3038 = t2992 + t3036;
  t3255 = -1.*t1779*t1503*t2956;
  t3260 = -1.*t1660*t3122;
  t3262 = t3255 + t3260;
  t2662 = -0.3455*t2398;
  t2705 = 0.35*t2422;
  t2719 = t2662 + t2705;
  t3268 = t2094*t3038;
  t3280 = t1049*t3262;
  t3289 = t3268 + t3280;
  t3322 = t1049*t3038;
  t3324 = -1.*t2094*t3262;
  t3330 = t3322 + t3324;
  t3482 = t1514*t2956;
  t3495 = t2996*t1707*t1904;
  t3497 = t3482 + t3495;
  t3453 = -1.*t2996*t1514*t1707;
  t3457 = t2956*t1904;
  t3458 = t3453 + t3457;
  t3520 = t1779*t2996*t1503;
  t3521 = -1.*t1660*t3497;
  t3527 = t3520 + t3521;
  t3539 = t2094*t3458;
  t3657 = t1049*t3527;
  t3660 = t3539 + t3657;
  t3726 = t1049*t3458;
  t3747 = -1.*t2094*t3527;
  t3749 = t3726 + t3747;
  p_output1[0]=-0.3455*t1466*t1503*t1514 + 0.072*t1660*t1707 + 0.072*t1503*t1842*t1904 - 0.19875*(t1660*t1707 - 1.*t1503*t1779*t1904) - 0.3455*t2094*t2221 + t2496*t2600 + t2719*t2748 - 0.3455*(-1.*t2422*t2600 + t2338*t2748) - 0.7*(t2338*t2600 + t2422*t2748) + var1[0];
  p_output1[1]=-0.072*t1503*t1660*t2956 - 0.3455*t1466*t3038 - 0.072*t1842*t3122 - 0.19875*(-1.*t1503*t1660*t2956 + t1779*t3122) - 0.3455*t2094*t3262 + t2496*t3289 + t2719*t3330 - 0.3455*(-1.*t2422*t3289 + t2338*t3330) - 0.7*(t2338*t3289 + t2422*t3330) + var1[1];
  p_output1[2]=0.072*t1503*t1660*t2996 - 0.3455*t1466*t3458 - 0.072*t1842*t3497 - 0.19875*(t1503*t1660*t2996 + t1779*t3497) - 0.3455*t2094*t3527 + t2496*t3660 + t2719*t3749 - 0.3455*(-1.*t2422*t3660 + t2338*t3749) - 0.7*(t2338*t3660 + t2422*t3749) + var1[2];
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

#include "RR_foot.hh"

namespace SymFunction
{

void RR_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
