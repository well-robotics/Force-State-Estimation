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
  double t1490;
  double t1620;
  double t1347;
  double t1523;
  double t1886;
  double t2060;
  double t2343;
  double t2374;
  double t2376;
  double t2502;
  double t2505;
  double t55;
  double t2754;
  double t2868;
  double t2884;
  double t1599;
  double t1935;
  double t1940;
  double t2899;
  double t3088;
  double t3125;
  double t3128;
  double t3214;
  double t3280;
  double t3293;
  double t3315;
  double t2992;
  double t3036;
  double t3038;
  double t289;
  double t818;
  double t2065;
  double t2122;
  double t3538;
  double t3539;
  double t3657;
  double t3828;
  double t3830;
  double t3835;
  double t3495;
  double t3502;
  double t3507;
  double t3147;
  double t3217;
  double t3246;
  double t3322;
  double t3341;
  double t3355;
  double t3870;
  double t3873;
  double t3879;
  double t3865;
  double t3867;
  double t3868;
  double t4114;
  double t4117;
  double t4124;
  double t4159;
  double t4179;
  double t4187;
  double t4137;
  double t4139;
  double t4150;
  double t4366;
  double t4378;
  double t4379;
  double t4494;
  double t4526;
  double t4529;
  double t4465;
  double t4468;
  double t4473;
  double t4754;
  double t4770;
  double t4773;
  double t4839;
  double t4840;
  double t4849;
  double t4782;
  double t4791;
  double t4828;
  double t5060;
  double t5061;
  double t5062;
  double t5047;
  double t5048;
  double t5054;
  double t5157;
  double t5173;
  double t5184;
  double t5320;
  double t5328;
  double t5343;
  double t5254;
  double t5278;
  double t5281;
  double t5453;
  double t5458;
  double t5468;
  double t5493;
  double t5504;
  double t5506;
  double t5481;
  double t5482;
  double t5483;
  double t5602;
  double t5632;
  double t5639;
  double t5756;
  double t5761;
  double t5767;
  double t2691;
  double t6009;
  double t6013;
  double t5552;
  double t5572;
  double t5578;
  double t6251;
  double t6256;
  double t6267;
  double t6298;
  double t6304;
  double t6318;
  double t3889;
  double t6440;
  double t6461;
  double t6463;
  double t5960;
  double t5970;
  double t5996;
  double t6562;
  double t6575;
  double t6577;
  double t6579;
  double t6582;
  double t6598;
  double t6603;
  double t6606;
  double t6608;
  double t6327;
  double t6700;
  double t6703;
  double t6704;
  double t6374;
  double t6706;
  double t6707;
  double t6710;
  double t6723;
  double t6725;
  double t3896;
  double t3901;
  double t6489;
  double t6611;
  double t6848;
  double t6856;
  double t6878;
  double t6657;
  t1490 = Cos(var1[5]);
  t1620 = Sin(var1[3]);
  t1347 = Cos(var1[3]);
  t1523 = Sin(var1[4]);
  t1886 = Sin(var1[5]);
  t2060 = Cos(var1[6]);
  t2343 = -1.*t1490*t1620;
  t2374 = -1.*t1347*t1523*t1886;
  t2376 = t2343 + t2374;
  t2502 = Cos(var1[4]);
  t2505 = Sin(var1[6]);
  t55 = Cos(var1[7]);
  t2754 = -1.*t1347*t2502*t2060;
  t2868 = -1.*t2376*t2505;
  t2884 = t2754 + t2868;
  t1599 = t1347*t1490*t1523;
  t1935 = -1.*t1620*t1886;
  t1940 = t1599 + t1935;
  t2899 = Sin(var1[7]);
  t3088 = Cos(var1[8]);
  t3125 = -1.*t3088;
  t3128 = 1. + t3125;
  t3214 = Sin(var1[8]);
  t3280 = t55*t1940;
  t3293 = -1.*t2884*t2899;
  t3315 = t3280 + t3293;
  t2992 = t55*t2884;
  t3036 = t1940*t2899;
  t3038 = t2992 + t3036;
  t289 = -1.*t55;
  t818 = 1. + t289;
  t2065 = -1.*t2060;
  t2122 = 1. + t2065;
  t3538 = t1347*t1490;
  t3539 = -1.*t1620*t1523*t1886;
  t3657 = t3538 + t3539;
  t3828 = -1.*t2502*t2060*t1620;
  t3830 = -1.*t3657*t2505;
  t3835 = t3828 + t3830;
  t3495 = t1490*t1620*t1523;
  t3502 = t1347*t1886;
  t3507 = t3495 + t3502;
  t3147 = -0.35*t3128;
  t3217 = 0.3455*t3214;
  t3246 = t3147 + t3217;
  t3322 = 0.3455*t3128;
  t3341 = 0.35*t3214;
  t3355 = t3322 + t3341;
  t3870 = t55*t3507;
  t3873 = -1.*t3835*t2899;
  t3879 = t3870 + t3873;
  t3865 = t55*t3835;
  t3867 = t3507*t2899;
  t3868 = t3865 + t3867;
  t4114 = t2502*t2060;
  t4117 = -1.*t1523*t1886*t2505;
  t4124 = t4114 + t4117;
  t4159 = -1.*t1490*t55*t1523;
  t4179 = -1.*t4124*t2899;
  t4187 = t4159 + t4179;
  t4137 = t55*t4124;
  t4139 = -1.*t1490*t1523*t2899;
  t4150 = t4137 + t4139;
  t4366 = t2060*t1620*t1523;
  t4378 = t2502*t1620*t1886*t2505;
  t4379 = t4366 + t4378;
  t4494 = t2502*t1490*t55*t1620;
  t4526 = -1.*t4379*t2899;
  t4529 = t4494 + t4526;
  t4465 = t55*t4379;
  t4468 = t2502*t1490*t1620*t2899;
  t4473 = t4465 + t4468;
  t4754 = -1.*t1347*t2060*t1523;
  t4770 = -1.*t1347*t2502*t1886*t2505;
  t4773 = t4754 + t4770;
  t4839 = -1.*t1347*t2502*t1490*t55;
  t4840 = -1.*t4773*t2899;
  t4849 = t4839 + t4840;
  t4782 = t55*t4773;
  t4791 = -1.*t1347*t2502*t1490*t2899;
  t4828 = t4782 + t4791;
  t5060 = -1.*t2502*t55*t1886;
  t5061 = -1.*t2502*t1490*t2505*t2899;
  t5062 = t5060 + t5061;
  t5047 = t2502*t1490*t55*t2505;
  t5048 = -1.*t2502*t1886*t2899;
  t5054 = t5047 + t5048;
  t5157 = -1.*t1490*t1620*t1523;
  t5173 = -1.*t1347*t1886;
  t5184 = t5157 + t5173;
  t5320 = t55*t3657;
  t5328 = t5184*t2505*t2899;
  t5343 = t5320 + t5328;
  t5254 = -1.*t55*t5184*t2505;
  t5278 = t3657*t2899;
  t5281 = t5254 + t5278;
  t5453 = t1490*t1620;
  t5458 = t1347*t1523*t1886;
  t5468 = t5453 + t5458;
  t5493 = t55*t5468;
  t5504 = t1940*t2505*t2899;
  t5506 = t5493 + t5504;
  t5481 = -1.*t55*t1940*t2505;
  t5482 = t5468*t2899;
  t5483 = t5481 + t5482;
  t5602 = t2502*t2060*t1886;
  t5632 = -1.*t1523*t2505;
  t5639 = t5602 + t5632;
  t5756 = -1.*t2060*t3657;
  t5761 = t2502*t1620*t2505;
  t5767 = t5756 + t5761;
  t2691 = -1.*t1347*t2502*t2505;
  t6009 = -1.*t2060*t5468;
  t6013 = t6009 + t2691;
  t5552 = t2060*t1523;
  t5572 = t2502*t1886*t2505;
  t5578 = t5552 + t5572;
  t6251 = t2502*t1490*t55;
  t6256 = -1.*t5578*t2899;
  t6267 = t6251 + t6256;
  t6298 = -1.*t55*t5578;
  t6304 = -1.*t2502*t1490*t2899;
  t6318 = t6298 + t6304;
  t3889 = t3088*t3879;
  t6440 = -1.*t55*t3835;
  t6461 = -1.*t3507*t2899;
  t6463 = t6440 + t6461;
  t5960 = t1347*t2502*t2060;
  t5970 = -1.*t5468*t2505;
  t5996 = t5960 + t5970;
  t6562 = -1.*t1347*t1490*t1523;
  t6575 = t1620*t1886;
  t6577 = t6562 + t6575;
  t6579 = t55*t6577;
  t6582 = -1.*t5996*t2899;
  t6598 = t6579 + t6582;
  t6603 = -1.*t55*t5996;
  t6606 = -1.*t6577*t2899;
  t6608 = t6603 + t6606;
  t6327 = t3088*t6267;
  t6700 = t55*t5578;
  t6703 = t2502*t1490*t2899;
  t6704 = t6700 + t6703;
  t6374 = -1.*t6267*t3214;
  t6706 = 0.3455*t3088;
  t6707 = -0.35*t3214;
  t6710 = t6706 + t6707;
  t6723 = 0.35*t3088;
  t6725 = t6723 + t3217;
  t3896 = -1.*t3868*t3214;
  t3901 = t3889 + t3896;
  t6489 = -1.*t3879*t3214;
  t6611 = t3088*t6598;
  t6848 = t55*t5996;
  t6856 = t6577*t2899;
  t6878 = t6848 + t6856;
  t6657 = -1.*t6598*t3214;
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
  p_output1[10]=-0.072*t2122*t2376 - 0.072*t1347*t2502*t2505 - 0.19875*(t2060*t2376 + t2691) + 0.3455*t2884*t2899 + t3038*t3246 + 0.3455*(-1.*t3038*t3214 + t3088*t3315) - 0.7*(t3038*t3088 + t3214*t3315) + t3315*t3355 + 0.3455*t1940*t818;
  p_output1[11]=-0.072*t1620*t2502*t2505 - 0.072*t2122*t3657 - 0.19875*(-1.*t1620*t2502*t2505 + t2060*t3657) + 0.3455*t2899*t3835 + t3246*t3868 + t3355*t3879 - 0.7*(t3088*t3868 + t3214*t3879) + 0.3455*t3901 + 0.3455*t3507*t818;
  p_output1[12]=-0.072*t1523*t1886*t2122 + 0.072*t2502*t2505 - 0.19875*(t1523*t1886*t2060 + t2502*t2505) + 0.3455*t2899*t4124 + t3246*t4150 + t3355*t4187 + 0.3455*(-1.*t3214*t4150 + t3088*t4187) - 0.7*(t3088*t4150 + t3214*t4187) - 0.3455*t1490*t1523*t818;
  p_output1[13]=0.072*t1620*t1886*t2122*t2502 + 0.072*t1523*t1620*t2505 - 0.19875*(-1.*t1620*t1886*t2060*t2502 + t1523*t1620*t2505) + 0.3455*t2899*t4379 + t3246*t4473 + t3355*t4529 + 0.3455*(-1.*t3214*t4473 + t3088*t4529) - 0.7*(t3088*t4473 + t3214*t4529) + 0.3455*t1490*t1620*t2502*t818;
  p_output1[14]=-0.072*t1347*t1886*t2122*t2502 - 0.072*t1347*t1523*t2505 - 0.19875*(t1347*t1886*t2060*t2502 - 1.*t1347*t1523*t2505) + 0.3455*t2899*t4773 + t3246*t4828 + t3355*t4849 + 0.3455*(-1.*t3214*t4828 + t3088*t4849) - 0.7*(t3088*t4828 + t3214*t4849) - 0.3455*t1347*t1490*t2502*t818;
  p_output1[15]=0.19875*t1490*t2060*t2502 + 0.072*t1490*t2122*t2502 + 0.3455*t1490*t2502*t2505*t2899 + t3246*t5054 + t3355*t5062 + 0.3455*(-1.*t3214*t5054 + t3088*t5062) - 0.7*(t3088*t5054 + t3214*t5062) - 0.3455*t1886*t2502*t818;
  p_output1[16]=-0.19875*t2060*t5184 - 0.072*t2122*t5184 - 0.3455*t2505*t2899*t5184 + t3246*t5281 + t3355*t5343 + 0.3455*(-1.*t3214*t5281 + t3088*t5343) - 0.7*(t3088*t5281 + t3214*t5343) + 0.3455*t3657*t818;
  p_output1[17]=-0.19875*t1940*t2060 - 0.072*t1940*t2122 - 0.3455*t1940*t2505*t2899 + t3246*t5483 + t3355*t5506 + 0.3455*(-1.*t3214*t5483 + t3088*t5506) - 0.7*(t3088*t5483 + t3214*t5506) + 0.3455*t5468*t818;
  p_output1[18]=0.072*t1523*t2060 + 0.072*t1886*t2502*t2505 - 0.19875*t5578 + 0.3455*t2899*t5639 - 1.*t2899*t3355*t5639 + t3246*t55*t5639 - 0.7*(-1.*t2899*t3214*t5639 + t3088*t55*t5639) + 0.3455*(-1.*t2899*t3088*t5639 - 1.*t3214*t55*t5639);
  p_output1[19]=-0.072*t1620*t2060*t2502 - 0.072*t2505*t3657 - 0.19875*t3835 + 0.3455*t2899*t5767 - 1.*t2899*t3355*t5767 + t3246*t55*t5767 - 0.7*(-1.*t2899*t3214*t5767 + t3088*t55*t5767) + 0.3455*(-1.*t2899*t3088*t5767 - 1.*t3214*t55*t5767);
  p_output1[20]=0.072*t1347*t2060*t2502 - 0.072*t2505*t5468 - 0.19875*t5996 + 0.3455*t2899*t6013 - 1.*t2899*t3355*t6013 + t3246*t55*t6013 - 0.7*(-1.*t2899*t3214*t6013 + t3088*t55*t6013) + 0.3455*(-1.*t2899*t3088*t6013 - 1.*t3214*t55*t6013);
  p_output1[21]=0.3455*t1490*t2502*t2899 + 0.3455*t55*t5578 + t3246*t6267 + t3355*t6318 - 0.7*(t3214*t6318 + t6327) + 0.3455*(t3088*t6318 + t6374);
  p_output1[22]=0.3455*t2899*t3507 + t3246*t3879 + 0.3455*t3835*t55 + t3355*t6463 - 0.7*(t3889 + t3214*t6463) + 0.3455*(t3088*t6463 + t6489);
  p_output1[23]=0.3455*t55*t5996 + 0.3455*t2899*t6577 + t3246*t6598 + t3355*t6608 - 0.7*(t3214*t6608 + t6611) + 0.3455*(t3088*t6608 + t6657);
  p_output1[24]=0.3455*(t6374 - 1.*t3088*t6704) - 0.7*(t6327 - 1.*t3214*t6704) + t6704*t6710 + t6267*t6725;
  p_output1[25]=-0.7*t3901 + 0.3455*(-1.*t3088*t3868 + t6489) + t3868*t6710 + t3879*t6725;
  p_output1[26]=t6598*t6725 + t6710*t6878 + 0.3455*(t6657 - 1.*t3088*t6878) - 0.7*(t6611 - 1.*t3214*t6878);
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
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=0;
  p_output1[48]=0;
  p_output1[49]=0;
  p_output1[50]=0;
  p_output1[51]=0;
  p_output1[52]=0;
  p_output1[53]=0;
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

#include "J_FR.hh"

namespace SymFunction
{

void J_FR_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
