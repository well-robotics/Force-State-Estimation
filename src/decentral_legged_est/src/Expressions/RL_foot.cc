/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:29 GMT-05:00
 */

#include "math2mat.hpp"
#include "mdefs.hpp"

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t621;
  double t1206;
  double t1223;
  double t1295;
  double t1352;
  double t948;
  double t1726;
  double t76;
  double t1772;
  double t1779;
  double t1810;
  double t1885;
  double t1894;
  double t1918;
  double t1974;
  double t2060;
  double t2063;
  double t2082;
  double t2237;
  double t2239;
  double t2310;
  double t231;
  double t353;
  double t2643;
  double t1309;
  double t1328;
  double t2678;
  double t2730;
  double t2738;
  double t2739;
  double t1940;
  double t2024;
  double t2026;
  double t2662;
  double t2705;
  double t2719;
  double t2831;
  double t2880;
  double t2890;
  double t2179;
  double t2200;
  double t2221;
  double t2899;
  double t2907;
  double t2914;
  double t2930;
  double t2936;
  double t2953;
  double t3252;
  double t3255;
  double t3259;
  double t3180;
  double t3217;
  double t3233;
  double t3289;
  double t3306;
  double t3308;
  double t3322;
  double t3324;
  double t3328;
  double t3341;
  double t3350;
  double t3354;
  t621 = Cos(var1[4]);
  t1206 = Sin(var1[15]);
  t1223 = Sin(var1[4]);
  t1295 = Cos(var1[15]);
  t1352 = Sin(var1[5]);
  t948 = Cos(var1[5]);
  t1726 = Sin(var1[16]);
  t76 = Cos(var1[16]);
  t1772 = t1295*t1223;
  t1779 = t621*t1206*t1352;
  t1810 = t1772 + t1779;
  t1885 = Cos(var1[17]);
  t1894 = -1.*t1885;
  t1918 = 1. + t1894;
  t1974 = Sin(var1[17]);
  t2060 = t621*t948*t1726;
  t2063 = t76*t1810;
  t2082 = t2060 + t2063;
  t2237 = t76*t621*t948;
  t2239 = -1.*t1726*t1810;
  t2310 = t2237 + t2239;
  t231 = -1.*t76;
  t353 = 1. + t231;
  t2643 = Sin(var1[3]);
  t1309 = -1.*t1295;
  t1328 = 1. + t1309;
  t2678 = Cos(var1[3]);
  t2730 = t2678*t948;
  t2738 = -1.*t2643*t1223*t1352;
  t2739 = t2730 + t2738;
  t1940 = -0.35*t1918;
  t2024 = -0.3455*t1974;
  t2026 = t1940 + t2024;
  t2662 = t948*t2643*t1223;
  t2705 = t2678*t1352;
  t2719 = t2662 + t2705;
  t2831 = -1.*t1295*t621*t2643;
  t2880 = -1.*t1206*t2739;
  t2890 = t2831 + t2880;
  t2179 = -0.3455*t1918;
  t2200 = 0.35*t1974;
  t2221 = t2179 + t2200;
  t2899 = t1726*t2719;
  t2907 = t76*t2890;
  t2914 = t2899 + t2907;
  t2930 = t76*t2719;
  t2936 = -1.*t1726*t2890;
  t2953 = t2930 + t2936;
  t3252 = t948*t2643;
  t3255 = t2678*t1223*t1352;
  t3259 = t3252 + t3255;
  t3180 = -1.*t2678*t948*t1223;
  t3217 = t2643*t1352;
  t3233 = t3180 + t3217;
  t3289 = t1295*t2678*t621;
  t3306 = -1.*t1206*t3259;
  t3308 = t3289 + t3306;
  t3322 = t1726*t3233;
  t3324 = t76*t3308;
  t3328 = t3322 + t3324;
  t3341 = t76*t3233;
  t3350 = -1.*t1726*t3308;
  t3354 = t3341 + t3350;
  p_output1[0]=-0.072*t1206*t1223 - 0.3455*t1726*t1810 + t2026*t2082 + t2221*t2310 - 0.3455*(-1.*t1974*t2082 + t1885*t2310) - 0.7*(t1885*t2082 + t1974*t2310) - 0.072*t1328*t1352*t621 + 0.19875*(t1206*t1223 - 1.*t1295*t1352*t621) - 0.3455*t353*t621*t948 + var1[0];
  p_output1[1]=0.072*t1328*t2739 - 0.3455*t1726*t2890 + t2026*t2914 + t2221*t2953 - 0.3455*(-1.*t1974*t2914 + t1885*t2953) - 0.7*(t1885*t2914 + t1974*t2953) - 0.3455*t2719*t353 + 0.072*t1206*t2643*t621 + 0.19875*(t1295*t2739 - 1.*t1206*t2643*t621) + var1[1];
  p_output1[2]=0.072*t1328*t3259 - 0.3455*t1726*t3308 + t2026*t3328 + t2221*t3354 - 0.3455*(-1.*t1974*t3328 + t1885*t3354) - 0.7*(t1885*t3328 + t1974*t3354) - 0.3455*t3233*t353 - 0.072*t1206*t2678*t621 + 0.19875*(t1295*t3259 + t1206*t2678*t621) + var1[2];
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

#include "RL_foot.hh"

namespace SymFunction
{

void RL_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
