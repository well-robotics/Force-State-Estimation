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
  double t187;
  double t647;
  double t912;
  double t1148;
  double t1174;
  double t481;
  double t1328;
  double t1352;
  double t1357;
  double t350;
  double t1469;
  double t1616;
  double t1617;
  double t1626;
  double t1672;
  double t1779;
  double t1789;
  double t1801;
  double t1514;
  double t1528;
  double t1558;
  double t544;
  double t577;
  double t711;
  double t719;
  double t2185;
  double t2140;
  double t2239;
  double t2274;
  double t2296;
  double t2432;
  double t2434;
  double t2436;
  double t2179;
  double t2200;
  double t2221;
  double t1660;
  double t1707;
  double t1726;
  double t1810;
  double t1832;
  double t1885;
  double t2600;
  double t2604;
  double t2644;
  double t2496;
  double t2537;
  double t2545;
  double t2880;
  double t2890;
  double t2891;
  double t2953;
  double t2954;
  double t2976;
  double t2767;
  double t2778;
  double t2825;
  double t3038;
  double t3040;
  double t3094;
  double t2992;
  double t2996;
  double t3025;
  t187 = Cos(var1[4]);
  t647 = Cos(var1[6]);
  t912 = Sin(var1[5]);
  t1148 = Sin(var1[4]);
  t1174 = Sin(var1[6]);
  t481 = Cos(var1[7]);
  t1328 = t647*t1148;
  t1352 = t187*t912*t1174;
  t1357 = t1328 + t1352;
  t350 = Cos(var1[5]);
  t1469 = Sin(var1[7]);
  t1616 = Cos(var1[8]);
  t1617 = -1.*t1616;
  t1626 = 1. + t1617;
  t1672 = Sin(var1[8]);
  t1779 = t187*t350*t481;
  t1789 = -1.*t1357*t1469;
  t1801 = t1779 + t1789;
  t1514 = t481*t1357;
  t1528 = t187*t350*t1469;
  t1558 = t1514 + t1528;
  t544 = -1.*t481;
  t577 = 1. + t544;
  t711 = -1.*t647;
  t719 = 1. + t711;
  t2185 = Cos(var1[3]);
  t2140 = Sin(var1[3]);
  t2239 = t2185*t350;
  t2274 = -1.*t2140*t1148*t912;
  t2296 = t2239 + t2274;
  t2432 = -1.*t187*t647*t2140;
  t2434 = -1.*t2296*t1174;
  t2436 = t2432 + t2434;
  t2179 = t350*t2140*t1148;
  t2200 = t2185*t912;
  t2221 = t2179 + t2200;
  t1660 = -0.35*t1626;
  t1707 = 0.3455*t1672;
  t1726 = t1660 + t1707;
  t1810 = 0.3455*t1626;
  t1832 = 0.35*t1672;
  t1885 = t1810 + t1832;
  t2600 = t481*t2221;
  t2604 = -1.*t2436*t1469;
  t2644 = t2600 + t2604;
  t2496 = t481*t2436;
  t2537 = t2221*t1469;
  t2545 = t2496 + t2537;
  t2880 = t350*t2140;
  t2890 = t2185*t1148*t912;
  t2891 = t2880 + t2890;
  t2953 = t2185*t187*t647;
  t2954 = -1.*t2891*t1174;
  t2976 = t2953 + t2954;
  t2767 = -1.*t2185*t350*t1148;
  t2778 = t2140*t912;
  t2825 = t2767 + t2778;
  t3038 = t481*t2825;
  t3040 = -1.*t2976*t1469;
  t3094 = t3038 + t3040;
  t2992 = t481*t2976;
  t2996 = t2825*t1469;
  t3025 = t2992 + t2996;
  p_output1[0]=0.072*t1148*t1174 + 0.3455*t1357*t1469 + t1558*t1726 + 0.3455*(-1.*t1558*t1672 + t1616*t1801) - 0.7*(t1558*t1616 + t1672*t1801) + t1801*t1885 + 0.3455*t187*t350*t577 + 0.072*t187*t719*t912 - 0.19875*(t1148*t1174 - 1.*t187*t647*t912) + var1[0];
  p_output1[1]=-0.072*t1174*t187*t2140 + 0.3455*t1469*t2436 + t1726*t2545 + t1885*t2644 + 0.3455*(-1.*t1672*t2545 + t1616*t2644) - 0.7*(t1616*t2545 + t1672*t2644) + 0.3455*t2221*t577 - 0.19875*(-1.*t1174*t187*t2140 + t2296*t647) - 0.072*t2296*t719 + var1[1];
  p_output1[2]=0.072*t1174*t187*t2185 + 0.3455*t1469*t2976 + t1726*t3025 + t1885*t3094 + 0.3455*(-1.*t1672*t3025 + t1616*t3094) - 0.7*(t1616*t3025 + t1672*t3094) + 0.3455*t2825*t577 - 0.19875*(t1174*t187*t2185 + t2891*t647) - 0.072*t2891*t719 + var1[2];
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

#include "FR_foot.hh"

namespace SymFunction
{

void FR_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
