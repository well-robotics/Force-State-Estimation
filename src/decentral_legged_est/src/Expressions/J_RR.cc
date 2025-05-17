/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:32 GMT-05:00
 */

#include "math2mat.hpp"
#include "mdefs.hpp"

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t25;
  double t623;
  double t894;
  double t757;
  double t908;
  double t92;
  double t99;
  double t1015;
  double t1243;
  double t1258;
  double t1322;
  double t1681;
  double t830;
  double t934;
  double t960;
  double t262;
  double t1684;
  double t1762;
  double t1870;
  double t1944;
  double t1982;
  double t2050;
  double t2097;
  double t2168;
  double t2199;
  double t2256;
  double t2719;
  double t2727;
  double t2728;
  double t389;
  double t619;
  double t1034;
  double t1155;
  double t3268;
  double t3280;
  double t3341;
  double t2052;
  double t2118;
  double t2160;
  double t3226;
  double t3230;
  double t3245;
  double t3414;
  double t3425;
  double t3427;
  double t2567;
  double t2573;
  double t2681;
  double t3495;
  double t3502;
  double t3508;
  double t3522;
  double t3531;
  double t3538;
  double t3867;
  double t3869;
  double t3870;
  double t3888;
  double t3889;
  double t3900;
  double t3957;
  double t4002;
  double t4008;
  double t4189;
  double t4197;
  double t4204;
  double t4221;
  double t4232;
  double t4243;
  double t4265;
  double t4269;
  double t4273;
  double t4551;
  double t4560;
  double t4563;
  double t4602;
  double t4633;
  double t4670;
  double t4700;
  double t4722;
  double t4728;
  double t5022;
  double t5024;
  double t5032;
  double t5038;
  double t5047;
  double t5056;
  double t5222;
  double t5247;
  double t5252;
  double t5362;
  double t5369;
  double t5390;
  double t5396;
  double t5397;
  double t5411;
  double t5511;
  double t5522;
  double t5535;
  double t5538;
  double t5540;
  double t5544;
  double t5548;
  double t5549;
  double t5550;
  double t5713;
  double t5720;
  double t5731;
  double t5945;
  double t5960;
  double t5984;
  double t1578;
  double t6223;
  double t6246;
  double t5777;
  double t5792;
  double t5810;
  double t6425;
  double t6447;
  double t6463;
  double t6481;
  double t6482;
  double t6524;
  double t6617;
  double t6622;
  double t6628;
  double t3630;
  double t6284;
  double t6293;
  double t6305;
  double t6731;
  double t6734;
  double t6753;
  double t6783;
  double t6790;
  double t6791;
  double t6817;
  double t6819;
  double t6826;
  double t6938;
  double t6939;
  double t6943;
  double t6537;
  double t6578;
  double t6933;
  double t6935;
  double t6936;
  double t6948;
  double t6951;
  double t3599;
  double t3694;
  double t6703;
  double t6993;
  double t6999;
  double t7002;
  double t6841;
  double t6898;
  t25 = Cos(var1[3]);
  t623 = Cos(var1[5]);
  t894 = Sin(var1[3]);
  t757 = Sin(var1[4]);
  t908 = Sin(var1[5]);
  t92 = Cos(var1[4]);
  t99 = Sin(var1[12]);
  t1015 = Cos(var1[12]);
  t1243 = -1.*t623*t894;
  t1258 = -1.*t25*t757*t908;
  t1322 = t1243 + t1258;
  t1681 = Sin(var1[13]);
  t830 = t25*t623*t757;
  t934 = -1.*t894*t908;
  t960 = t830 + t934;
  t262 = Cos(var1[13]);
  t1684 = -1.*t1015*t25*t92;
  t1762 = -1.*t99*t1322;
  t1870 = t1684 + t1762;
  t1944 = Cos(var1[14]);
  t1982 = -1.*t1944;
  t2050 = 1. + t1982;
  t2097 = Sin(var1[14]);
  t2168 = t1681*t960;
  t2199 = t262*t1870;
  t2256 = t2168 + t2199;
  t2719 = t262*t960;
  t2727 = -1.*t1681*t1870;
  t2728 = t2719 + t2727;
  t389 = -1.*t262;
  t619 = 1. + t389;
  t1034 = -1.*t1015;
  t1155 = 1. + t1034;
  t3268 = t25*t623;
  t3280 = -1.*t894*t757*t908;
  t3341 = t3268 + t3280;
  t2052 = -0.35*t2050;
  t2118 = -0.3455*t2097;
  t2160 = t2052 + t2118;
  t3226 = t623*t894*t757;
  t3230 = t25*t908;
  t3245 = t3226 + t3230;
  t3414 = -1.*t1015*t92*t894;
  t3425 = -1.*t99*t3341;
  t3427 = t3414 + t3425;
  t2567 = -0.3455*t2050;
  t2573 = 0.35*t2097;
  t2681 = t2567 + t2573;
  t3495 = t1681*t3245;
  t3502 = t262*t3427;
  t3508 = t3495 + t3502;
  t3522 = t262*t3245;
  t3531 = -1.*t1681*t3427;
  t3538 = t3522 + t3531;
  t3867 = t1015*t92;
  t3869 = -1.*t99*t757*t908;
  t3870 = t3867 + t3869;
  t3888 = -1.*t623*t1681*t757;
  t3889 = t262*t3870;
  t3900 = t3888 + t3889;
  t3957 = -1.*t262*t623*t757;
  t4002 = -1.*t1681*t3870;
  t4008 = t3957 + t4002;
  t4189 = t1015*t894*t757;
  t4197 = t92*t99*t894*t908;
  t4204 = t4189 + t4197;
  t4221 = t92*t623*t1681*t894;
  t4232 = t262*t4204;
  t4243 = t4221 + t4232;
  t4265 = t262*t92*t623*t894;
  t4269 = -1.*t1681*t4204;
  t4273 = t4265 + t4269;
  t4551 = -1.*t1015*t25*t757;
  t4560 = -1.*t25*t92*t99*t908;
  t4563 = t4551 + t4560;
  t4602 = -1.*t25*t92*t623*t1681;
  t4633 = t262*t4563;
  t4670 = t4602 + t4633;
  t4700 = -1.*t262*t25*t92*t623;
  t4722 = -1.*t1681*t4563;
  t4728 = t4700 + t4722;
  t5022 = -1.*t92*t623*t99*t1681;
  t5024 = -1.*t262*t92*t908;
  t5032 = t5022 + t5024;
  t5038 = t262*t92*t623*t99;
  t5047 = -1.*t92*t1681*t908;
  t5056 = t5038 + t5047;
  t5222 = -1.*t623*t894*t757;
  t5247 = -1.*t25*t908;
  t5252 = t5222 + t5247;
  t5362 = t99*t1681*t5252;
  t5369 = t262*t3341;
  t5390 = t5362 + t5369;
  t5396 = -1.*t262*t99*t5252;
  t5397 = t1681*t3341;
  t5411 = t5396 + t5397;
  t5511 = t623*t894;
  t5522 = t25*t757*t908;
  t5535 = t5511 + t5522;
  t5538 = t99*t1681*t960;
  t5540 = t262*t5535;
  t5544 = t5538 + t5540;
  t5548 = -1.*t262*t99*t960;
  t5549 = t1681*t5535;
  t5550 = t5548 + t5549;
  t5713 = -1.*t99*t757;
  t5720 = t1015*t92*t908;
  t5731 = t5713 + t5720;
  t5945 = t92*t99*t894;
  t5960 = -1.*t1015*t3341;
  t5984 = t5945 + t5960;
  t1578 = -1.*t25*t92*t99;
  t6223 = -1.*t1015*t5535;
  t6246 = t1578 + t6223;
  t5777 = t1015*t757;
  t5792 = t92*t99*t908;
  t5810 = t5777 + t5792;
  t6425 = -1.*t92*t623*t1681;
  t6447 = -1.*t262*t5810;
  t6463 = t6425 + t6447;
  t6481 = t262*t92*t623;
  t6482 = -1.*t1681*t5810;
  t6524 = t6481 + t6482;
  t6617 = -1.*t1681*t3245;
  t6622 = -1.*t262*t3427;
  t6628 = t6617 + t6622;
  t3630 = t1944*t3538;
  t6284 = t1015*t25*t92;
  t6293 = -1.*t99*t5535;
  t6305 = t6284 + t6293;
  t6731 = -1.*t25*t623*t757;
  t6734 = t894*t908;
  t6753 = t6731 + t6734;
  t6783 = -1.*t1681*t6753;
  t6790 = -1.*t262*t6305;
  t6791 = t6783 + t6790;
  t6817 = t262*t6753;
  t6819 = -1.*t1681*t6305;
  t6826 = t6817 + t6819;
  t6938 = t92*t623*t1681;
  t6939 = t262*t5810;
  t6943 = t6938 + t6939;
  t6537 = t1944*t6524;
  t6578 = -1.*t2097*t6524;
  t6933 = -0.3455*t1944;
  t6935 = -0.35*t2097;
  t6936 = t6933 + t6935;
  t6948 = 0.35*t1944;
  t6951 = t6948 + t2118;
  t3599 = -1.*t2097*t3508;
  t3694 = t3599 + t3630;
  t6703 = -1.*t2097*t3538;
  t6993 = t1681*t6753;
  t6999 = t262*t6305;
  t7002 = t6993 + t6999;
  t6841 = t1944*t6826;
  t6898 = -1.*t2097*t6826;
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
  p_output1[10]=-0.072*t1155*t1322 - 0.19875*(t1015*t1322 + t1578) - 0.3455*t1681*t1870 + t2160*t2256 + t2681*t2728 - 0.3455*(-1.*t2097*t2256 + t1944*t2728) - 0.7*(t1944*t2256 + t2097*t2728) - 0.3455*t619*t960 - 0.072*t25*t92*t99;
  p_output1[11]=-0.072*t1155*t3341 - 0.3455*t1681*t3427 + t2160*t3508 + t2681*t3538 - 0.7*(t1944*t3508 + t2097*t3538) - 0.3455*t3694 - 0.3455*t3245*t619 - 0.072*t894*t92*t99 - 0.19875*(t1015*t3341 - 1.*t894*t92*t99);
  p_output1[12]=-0.3455*t1681*t3870 + t2160*t3900 + t2681*t4008 - 0.3455*(-1.*t2097*t3900 + t1944*t4008) - 0.7*(t1944*t3900 + t2097*t4008) + 0.3455*t619*t623*t757 - 0.072*t1155*t757*t908 + 0.072*t92*t99 - 0.19875*(t1015*t757*t908 + t92*t99);
  p_output1[13]=-0.3455*t1681*t4204 + t2160*t4243 + t2681*t4273 - 0.3455*(-1.*t2097*t4243 + t1944*t4273) - 0.7*(t1944*t4243 + t2097*t4273) - 0.3455*t619*t623*t894*t92 + 0.072*t1155*t894*t908*t92 + 0.072*t757*t894*t99 - 0.19875*(-1.*t1015*t894*t908*t92 + t757*t894*t99);
  p_output1[14]=-0.3455*t1681*t4563 + t2160*t4670 + t2681*t4728 - 0.3455*(-1.*t2097*t4670 + t1944*t4728) - 0.7*(t1944*t4670 + t2097*t4728) + 0.3455*t25*t619*t623*t92 - 0.072*t1155*t25*t908*t92 - 0.072*t25*t757*t99 - 0.19875*(t1015*t25*t908*t92 - 1.*t25*t757*t99);
  p_output1[15]=t2681*t5032 + t2160*t5056 - 0.7*(t2097*t5032 + t1944*t5056) - 0.3455*(t1944*t5032 - 1.*t2097*t5056) + 0.19875*t1015*t623*t92 + 0.072*t1155*t623*t92 + 0.3455*t619*t908*t92 - 0.3455*t1681*t623*t92*t99;
  p_output1[16]=-0.19875*t1015*t5252 - 0.072*t1155*t5252 + t2681*t5390 + t2160*t5411 - 0.7*(t2097*t5390 + t1944*t5411) - 0.3455*(t1944*t5390 - 1.*t2097*t5411) - 0.3455*t3341*t619 + 0.3455*t1681*t5252*t99;
  p_output1[17]=t2681*t5544 + t2160*t5550 - 0.7*(t2097*t5544 + t1944*t5550) - 0.3455*(t1944*t5544 - 1.*t2097*t5550) - 0.3455*t5535*t619 - 0.19875*t1015*t960 - 0.072*t1155*t960 + 0.3455*t1681*t960*t99;
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
  p_output1[36]=-0.3455*t1681*t5731 + t2160*t262*t5731 - 1.*t1681*t2681*t5731 - 0.7*(-1.*t1681*t2097*t5731 + t1944*t262*t5731) - 0.3455*(-1.*t1681*t1944*t5731 - 1.*t2097*t262*t5731) - 0.19875*t5810 + 0.072*t1015*t757 + 0.072*t908*t92*t99;
  p_output1[37]=-0.19875*t3427 - 0.3455*t1681*t5984 + t2160*t262*t5984 - 1.*t1681*t2681*t5984 - 0.7*(-1.*t1681*t2097*t5984 + t1944*t262*t5984) - 0.3455*(-1.*t1681*t1944*t5984 - 1.*t2097*t262*t5984) - 0.072*t1015*t894*t92 - 0.072*t3341*t99;
  p_output1[38]=-0.3455*t1681*t6246 + t2160*t262*t6246 - 1.*t1681*t2681*t6246 - 0.7*(-1.*t1681*t2097*t6246 + t1944*t262*t6246) - 0.3455*(-1.*t1681*t1944*t6246 - 1.*t2097*t262*t6246) - 0.19875*t6305 + 0.072*t1015*t25*t92 - 0.072*t5535*t99;
  p_output1[39]=-0.3455*t262*t5810 + t2681*t6463 + t2160*t6524 - 0.7*(t2097*t6463 + t6537) - 0.3455*(t1944*t6463 + t6578) - 0.3455*t1681*t623*t92;
  p_output1[40]=-0.3455*t1681*t3245 - 0.3455*t262*t3427 + t2160*t3538 + t2681*t6628 - 0.7*(t3630 + t2097*t6628) - 0.3455*(t1944*t6628 + t6703);
  p_output1[41]=-0.3455*t262*t6305 - 0.3455*t1681*t6753 + t2681*t6791 + t2160*t6826 - 0.7*(t2097*t6791 + t6841) - 0.3455*(t1944*t6791 + t6898);
  p_output1[42]=t6936*t6943 - 0.3455*(t6578 - 1.*t1944*t6943) - 0.7*(t6537 - 1.*t2097*t6943) + t6524*t6951;
  p_output1[43]=-0.7*t3694 - 0.3455*(-1.*t1944*t3508 + t6703) + t3508*t6936 + t3538*t6951;
  p_output1[44]=t6826*t6951 + t6936*t7002 - 0.3455*(t6898 - 1.*t1944*t7002) - 0.7*(t6841 - 1.*t2097*t7002);
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

#include "J_RR.hh"

namespace SymFunction
{

void J_RR_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
