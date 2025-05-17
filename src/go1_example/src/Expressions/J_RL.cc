/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:31 GMT-05:00
 */

#include "math2mat.hpp"
#include "mdefs.hpp"

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t167;
  double t805;
  double t1185;
  double t842;
  double t1215;
  double t194;
  double t374;
  double t1247;
  double t1366;
  double t1385;
  double t1520;
  double t1690;
  double t908;
  double t1225;
  double t1233;
  double t486;
  double t1714;
  double t1740;
  double t1750;
  double t1828;
  double t1857;
  double t1916;
  double t1953;
  double t2168;
  double t2199;
  double t2410;
  double t2681;
  double t2685;
  double t2718;
  double t549;
  double t705;
  double t1255;
  double t1270;
  double t3359;
  double t3400;
  double t3404;
  double t1944;
  double t1982;
  double t2052;
  double t3268;
  double t3280;
  double t3322;
  double t3452;
  double t3457;
  double t3482;
  double t2567;
  double t2573;
  double t2677;
  double t3502;
  double t3510;
  double t3536;
  double t3599;
  double t3651;
  double t3690;
  double t3900;
  double t3901;
  double t3906;
  double t3946;
  double t3957;
  double t3970;
  double t4018;
  double t4029;
  double t4043;
  double t4221;
  double t4265;
  double t4266;
  double t4273;
  double t4288;
  double t4289;
  double t4304;
  double t4306;
  double t4323;
  double t4585;
  double t4592;
  double t4599;
  double t4608;
  double t4633;
  double t4665;
  double t4675;
  double t4678;
  double t4683;
  double t4963;
  double t4976;
  double t4977;
  double t4999;
  double t5004;
  double t5008;
  double t5133;
  double t5134;
  double t5136;
  double t5247;
  double t5252;
  double t5278;
  double t5320;
  double t5328;
  double t5359;
  double t5469;
  double t5470;
  double t5473;
  double t5486;
  double t5493;
  double t5504;
  double t5515;
  double t5522;
  double t5529;
  double t5642;
  double t5653;
  double t5660;
  double t5872;
  double t5873;
  double t5874;
  double t1599;
  double t6076;
  double t6082;
  double t5686;
  double t5688;
  double t5689;
  double t6393;
  double t6398;
  double t6413;
  double t6465;
  double t6466;
  double t6468;
  double t6601;
  double t6608;
  double t6610;
  double t3737;
  double t6223;
  double t6225;
  double t6228;
  double t6700;
  double t6703;
  double t6706;
  double t6718;
  double t6723;
  double t6725;
  double t6734;
  double t6742;
  double t6750;
  double t6841;
  double t6848;
  double t6856;
  double t6493;
  double t6546;
  double t6826;
  double t6830;
  double t6834;
  double t6898;
  double t6902;
  double t3726;
  double t3750;
  double t6664;
  double t6966;
  double t6970;
  double t6982;
  double t6778;
  double t6794;
  t167 = Cos(var1[3]);
  t805 = Cos(var1[5]);
  t1185 = Sin(var1[3]);
  t842 = Sin(var1[4]);
  t1215 = Sin(var1[5]);
  t194 = Cos(var1[4]);
  t374 = Sin(var1[15]);
  t1247 = Cos(var1[15]);
  t1366 = -1.*t805*t1185;
  t1385 = -1.*t167*t842*t1215;
  t1520 = t1366 + t1385;
  t1690 = Sin(var1[16]);
  t908 = t167*t805*t842;
  t1225 = -1.*t1185*t1215;
  t1233 = t908 + t1225;
  t486 = Cos(var1[16]);
  t1714 = -1.*t1247*t167*t194;
  t1740 = -1.*t374*t1520;
  t1750 = t1714 + t1740;
  t1828 = Cos(var1[17]);
  t1857 = -1.*t1828;
  t1916 = 1. + t1857;
  t1953 = Sin(var1[17]);
  t2168 = t1690*t1233;
  t2199 = t486*t1750;
  t2410 = t2168 + t2199;
  t2681 = t486*t1233;
  t2685 = -1.*t1690*t1750;
  t2718 = t2681 + t2685;
  t549 = -1.*t486;
  t705 = 1. + t549;
  t1255 = -1.*t1247;
  t1270 = 1. + t1255;
  t3359 = t167*t805;
  t3400 = -1.*t1185*t842*t1215;
  t3404 = t3359 + t3400;
  t1944 = -0.35*t1916;
  t1982 = -0.3455*t1953;
  t2052 = t1944 + t1982;
  t3268 = t805*t1185*t842;
  t3280 = t167*t1215;
  t3322 = t3268 + t3280;
  t3452 = -1.*t1247*t194*t1185;
  t3457 = -1.*t374*t3404;
  t3482 = t3452 + t3457;
  t2567 = -0.3455*t1916;
  t2573 = 0.35*t1953;
  t2677 = t2567 + t2573;
  t3502 = t1690*t3322;
  t3510 = t486*t3482;
  t3536 = t3502 + t3510;
  t3599 = t486*t3322;
  t3651 = -1.*t1690*t3482;
  t3690 = t3599 + t3651;
  t3900 = t1247*t194;
  t3901 = -1.*t374*t842*t1215;
  t3906 = t3900 + t3901;
  t3946 = -1.*t805*t1690*t842;
  t3957 = t486*t3906;
  t3970 = t3946 + t3957;
  t4018 = -1.*t486*t805*t842;
  t4029 = -1.*t1690*t3906;
  t4043 = t4018 + t4029;
  t4221 = t1247*t1185*t842;
  t4265 = t194*t374*t1185*t1215;
  t4266 = t4221 + t4265;
  t4273 = t194*t805*t1690*t1185;
  t4288 = t486*t4266;
  t4289 = t4273 + t4288;
  t4304 = t486*t194*t805*t1185;
  t4306 = -1.*t1690*t4266;
  t4323 = t4304 + t4306;
  t4585 = -1.*t1247*t167*t842;
  t4592 = -1.*t167*t194*t374*t1215;
  t4599 = t4585 + t4592;
  t4608 = -1.*t167*t194*t805*t1690;
  t4633 = t486*t4599;
  t4665 = t4608 + t4633;
  t4675 = -1.*t486*t167*t194*t805;
  t4678 = -1.*t1690*t4599;
  t4683 = t4675 + t4678;
  t4963 = -1.*t194*t805*t374*t1690;
  t4976 = -1.*t486*t194*t1215;
  t4977 = t4963 + t4976;
  t4999 = t486*t194*t805*t374;
  t5004 = -1.*t194*t1690*t1215;
  t5008 = t4999 + t5004;
  t5133 = -1.*t805*t1185*t842;
  t5134 = -1.*t167*t1215;
  t5136 = t5133 + t5134;
  t5247 = t374*t1690*t5136;
  t5252 = t486*t3404;
  t5278 = t5247 + t5252;
  t5320 = -1.*t486*t374*t5136;
  t5328 = t1690*t3404;
  t5359 = t5320 + t5328;
  t5469 = t805*t1185;
  t5470 = t167*t842*t1215;
  t5473 = t5469 + t5470;
  t5486 = t374*t1690*t1233;
  t5493 = t486*t5473;
  t5504 = t5486 + t5493;
  t5515 = -1.*t486*t374*t1233;
  t5522 = t1690*t5473;
  t5529 = t5515 + t5522;
  t5642 = -1.*t374*t842;
  t5653 = t1247*t194*t1215;
  t5660 = t5642 + t5653;
  t5872 = t194*t374*t1185;
  t5873 = -1.*t1247*t3404;
  t5874 = t5872 + t5873;
  t1599 = -1.*t167*t194*t374;
  t6076 = -1.*t1247*t5473;
  t6082 = t1599 + t6076;
  t5686 = t1247*t842;
  t5688 = t194*t374*t1215;
  t5689 = t5686 + t5688;
  t6393 = -1.*t194*t805*t1690;
  t6398 = -1.*t486*t5689;
  t6413 = t6393 + t6398;
  t6465 = t486*t194*t805;
  t6466 = -1.*t1690*t5689;
  t6468 = t6465 + t6466;
  t6601 = -1.*t1690*t3322;
  t6608 = -1.*t486*t3482;
  t6610 = t6601 + t6608;
  t3737 = t1828*t3690;
  t6223 = t1247*t167*t194;
  t6225 = -1.*t374*t5473;
  t6228 = t6223 + t6225;
  t6700 = -1.*t167*t805*t842;
  t6703 = t1185*t1215;
  t6706 = t6700 + t6703;
  t6718 = -1.*t1690*t6706;
  t6723 = -1.*t486*t6228;
  t6725 = t6718 + t6723;
  t6734 = t486*t6706;
  t6742 = -1.*t1690*t6228;
  t6750 = t6734 + t6742;
  t6841 = t194*t805*t1690;
  t6848 = t486*t5689;
  t6856 = t6841 + t6848;
  t6493 = t1828*t6468;
  t6546 = -1.*t1953*t6468;
  t6826 = -0.3455*t1828;
  t6830 = -0.35*t1953;
  t6834 = t6826 + t6830;
  t6898 = 0.35*t1828;
  t6902 = t6898 + t1982;
  t3726 = -1.*t1953*t3536;
  t3750 = t3726 + t3737;
  t6664 = -1.*t1953*t3690;
  t6966 = t1690*t6706;
  t6970 = t486*t6228;
  t6982 = t6966 + t6970;
  t6778 = t1828*t6750;
  t6794 = -1.*t1953*t6750;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=1.;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=1.;
  p_output1[9]=0;
  p_output1[10]=0.072*t1270*t1520 + 0.19875*(t1247*t1520 + t1599) - 0.3455*t1690*t1750 + t2052*t2410 + t2677*t2718 - 0.3455*(-1.*t1953*t2410 + t1828*t2718) - 0.7*(t1828*t2410 + t1953*t2718) + 0.072*t167*t194*t374 - 0.3455*t1233*t705;
  p_output1[11]=0.072*t1270*t3404 - 0.3455*t1690*t3482 + t2052*t3536 + t2677*t3690 - 0.7*(t1828*t3536 + t1953*t3690) + 0.072*t1185*t194*t374 + 0.19875*(t1247*t3404 - 1.*t1185*t194*t374) - 0.3455*t3750 - 0.3455*t3322*t705;
  p_output1[12]=-0.072*t194*t374 - 0.3455*t1690*t3906 + t2052*t3970 + t2677*t4043 - 0.3455*(-1.*t1953*t3970 + t1828*t4043) - 0.7*(t1828*t3970 + t1953*t4043) + 0.072*t1215*t1270*t842 + 0.3455*t705*t805*t842 + 0.19875*(t194*t374 + t1215*t1247*t842);
  p_output1[13]=-0.072*t1185*t1215*t1270*t194 - 0.3455*t1690*t4266 + t2052*t4289 + t2677*t4323 - 0.3455*(-1.*t1953*t4289 + t1828*t4323) - 0.7*(t1828*t4289 + t1953*t4323) - 0.3455*t1185*t194*t705*t805 - 0.072*t1185*t374*t842 + 0.19875*(-1.*t1185*t1215*t1247*t194 + t1185*t374*t842);
  p_output1[14]=0.072*t1215*t1270*t167*t194 - 0.3455*t1690*t4599 + t2052*t4665 + t2677*t4683 - 0.3455*(-1.*t1953*t4665 + t1828*t4683) - 0.7*(t1828*t4665 + t1953*t4683) + 0.3455*t167*t194*t705*t805 + 0.072*t167*t374*t842 + 0.19875*(t1215*t1247*t167*t194 - 1.*t167*t374*t842);
  p_output1[15]=t2677*t4977 + t2052*t5008 - 0.7*(t1953*t4977 + t1828*t5008) - 0.3455*(t1828*t4977 - 1.*t1953*t5008) + 0.3455*t1215*t194*t705 - 0.19875*t1247*t194*t805 - 0.072*t1270*t194*t805 - 0.3455*t1690*t194*t374*t805;
  p_output1[16]=0.19875*t1247*t5136 + 0.072*t1270*t5136 + 0.3455*t1690*t374*t5136 + t2677*t5278 + t2052*t5359 - 0.7*(t1953*t5278 + t1828*t5359) - 0.3455*(t1828*t5278 - 1.*t1953*t5359) - 0.3455*t3404*t705;
  p_output1[17]=0.19875*t1233*t1247 + 0.072*t1233*t1270 + 0.3455*t1233*t1690*t374 + t2677*t5504 + t2052*t5529 - 0.7*(t1953*t5504 + t1828*t5529) - 0.3455*(t1828*t5504 - 1.*t1953*t5529) - 0.3455*t5473*t705;
  p_output1[18]=0;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=0;
  p_output1[22]=0;
  p_output1[23]=0;
  p_output1[24]=0;
  p_output1[25]=0;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=0;
  p_output1[29]=0;
  p_output1[30]=0;
  p_output1[31]=0;
  p_output1[32]=0;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=0;
  p_output1[36]=0;
  p_output1[37]=0;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0;
  p_output1[41]=0;
  p_output1[42]=0;
  p_output1[43]=0;
  p_output1[44]=0;
  p_output1[45]=-0.072*t1215*t194*t374 - 0.3455*t1690*t5660 - 1.*t1690*t2677*t5660 + t2052*t486*t5660 - 0.7*(-1.*t1690*t1953*t5660 + t1828*t486*t5660) - 0.3455*(-1.*t1690*t1828*t5660 - 1.*t1953*t486*t5660) + 0.19875*t5689 - 0.072*t1247*t842;
  p_output1[46]=0.072*t1185*t1247*t194 + 0.19875*t3482 + 0.072*t3404*t374 - 0.3455*t1690*t5874 - 1.*t1690*t2677*t5874 + t2052*t486*t5874 - 0.7*(-1.*t1690*t1953*t5874 + t1828*t486*t5874) - 0.3455*(-1.*t1690*t1828*t5874 - 1.*t1953*t486*t5874);
  p_output1[47]=-0.072*t1247*t167*t194 + 0.072*t374*t5473 - 0.3455*t1690*t6082 - 1.*t1690*t2677*t6082 + t2052*t486*t6082 - 0.7*(-1.*t1690*t1953*t6082 + t1828*t486*t6082) - 0.3455*(-1.*t1690*t1828*t6082 - 1.*t1953*t486*t6082) + 0.19875*t6228;
  p_output1[48]=-0.3455*t486*t5689 + t2677*t6413 + t2052*t6468 - 0.7*(t1953*t6413 + t6493) - 0.3455*(t1828*t6413 + t6546) - 0.3455*t1690*t194*t805;
  p_output1[49]=-0.3455*t1690*t3322 + t2052*t3690 - 0.3455*t3482*t486 + t2677*t6610 - 0.7*(t3737 + t1953*t6610) - 0.3455*(t1828*t6610 + t6664);
  p_output1[50]=-0.3455*t486*t6228 - 0.3455*t1690*t6706 + t2677*t6725 + t2052*t6750 - 0.7*(t1953*t6725 + t6778) - 0.3455*(t1828*t6725 + t6794);
  p_output1[51]=t6834*t6856 - 0.3455*(t6546 - 1.*t1828*t6856) - 0.7*(t6493 - 1.*t1953*t6856) + t6468*t6902;
  p_output1[52]=-0.7*t3750 - 0.3455*(-1.*t1828*t3536 + t6664) + t3536*t6834 + t3690*t6902;
  p_output1[53]=t6750*t6902 + t6834*t6982 - 0.3455*(t6794 - 1.*t1828*t6982) - 0.7*(t6778 - 1.*t1953*t6982);
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 18, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "J_RL.hh"

namespace SymFunction
{

void J_RL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
