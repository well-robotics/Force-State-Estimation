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
  double t932;
  double t2030;
  double t2063;
  double t2041;
  double t2179;
  double t1528;
  double t1772;
  double t2478;
  double t2578;
  double t2653;
  double t2662;
  double t2765;
  double t2060;
  double t2200;
  double t2237;
  double t2021;
  double t2778;
  double t2825;
  double t2897;
  double t2907;
  double t2913;
  double t2924;
  double t2982;
  double t3038;
  double t3088;
  double t3110;
  double t3246;
  double t3252;
  double t3263;
  double t2024;
  double t2026;
  double t2496;
  double t2537;
  double t3457;
  double t3458;
  double t3481;
  double t2930;
  double t2992;
  double t3036;
  double t3425;
  double t3433;
  double t3442;
  double t3512;
  double t3514;
  double t3520;
  double t3214;
  double t3217;
  double t3233;
  double t3539;
  double t3657;
  double t3686;
  double t3759;
  double t3763;
  double t3768;
  double t3901;
  double t3908;
  double t3932;
  double t3935;
  double t3957;
  double t3959;
  double t3978;
  double t3982;
  double t4013;
  double t4159;
  double t4179;
  double t4181;
  double t4189;
  double t4197;
  double t4203;
  double t4221;
  double t4226;
  double t4234;
  double t4390;
  double t4414;
  double t4458;
  double t4467;
  double t4468;
  double t4469;
  double t4480;
  double t4481;
  double t4492;
  double t4687;
  double t4696;
  double t4699;
  double t4728;
  double t4734;
  double t4750;
  double t4869;
  double t4885;
  double t4887;
  double t4934;
  double t4984;
  double t4997;
  double t5002;
  double t5004;
  double t5034;
  double t5134;
  double t5136;
  double t5142;
  double t5151;
  double t5157;
  double t5169;
  double t5182;
  double t5184;
  double t5195;
  double t5363;
  double t5369;
  double t5376;
  double t5493;
  double t5496;
  double t5502;
  double t2705;
  double t5555;
  double t5569;
  double t5411;
  double t5417;
  double t5418;
  double t5723;
  double t5727;
  double t5731;
  double t5761;
  double t5770;
  double t5805;
  double t5900;
  double t5910;
  double t5925;
  double t3824;
  double t5602;
  double t5622;
  double t5629;
  double t6032;
  double t6034;
  double t6040;
  double t6073;
  double t6079;
  double t6082;
  double t6157;
  double t6162;
  double t6168;
  double t6331;
  double t6333;
  double t6344;
  double t5821;
  double t5870;
  double t6320;
  double t6327;
  double t6330;
  double t6374;
  double t6376;
  double t3814;
  double t3828;
  double t6015;
  double t6537;
  double t6546;
  double t6548;
  double t6252;
  double t6297;
  t932 = Cos(var1[3]);
  t2030 = Cos(var1[5]);
  t2063 = Sin(var1[3]);
  t2041 = Sin(var1[4]);
  t2179 = Sin(var1[5]);
  t1528 = Cos(var1[4]);
  t1772 = Sin(var1[9]);
  t2478 = Cos(var1[9]);
  t2578 = -1.*t2030*t2063;
  t2653 = -1.*t932*t2041*t2179;
  t2662 = t2578 + t2653;
  t2765 = Sin(var1[10]);
  t2060 = t932*t2030*t2041;
  t2200 = -1.*t2063*t2179;
  t2237 = t2060 + t2200;
  t2021 = Cos(var1[10]);
  t2778 = -1.*t2478*t932*t1528;
  t2825 = -1.*t1772*t2662;
  t2897 = t2778 + t2825;
  t2907 = Cos(var1[11]);
  t2913 = -1.*t2907;
  t2924 = 1. + t2913;
  t2982 = Sin(var1[11]);
  t3038 = t2765*t2237;
  t3088 = t2021*t2897;
  t3110 = t3038 + t3088;
  t3246 = t2021*t2237;
  t3252 = -1.*t2765*t2897;
  t3263 = t3246 + t3252;
  t2024 = -1.*t2021;
  t2026 = 1. + t2024;
  t2496 = -1.*t2478;
  t2537 = 1. + t2496;
  t3457 = t932*t2030;
  t3458 = -1.*t2063*t2041*t2179;
  t3481 = t3457 + t3458;
  t2930 = -0.35*t2924;
  t2992 = 0.3455*t2982;
  t3036 = t2930 + t2992;
  t3425 = t2030*t2063*t2041;
  t3433 = t932*t2179;
  t3442 = t3425 + t3433;
  t3512 = -1.*t2478*t1528*t2063;
  t3514 = -1.*t1772*t3481;
  t3520 = t3512 + t3514;
  t3214 = 0.3455*t2924;
  t3217 = 0.35*t2982;
  t3233 = t3214 + t3217;
  t3539 = t2765*t3442;
  t3657 = t2021*t3520;
  t3686 = t3539 + t3657;
  t3759 = t2021*t3442;
  t3763 = -1.*t2765*t3520;
  t3768 = t3759 + t3763;
  t3901 = t2478*t1528;
  t3908 = -1.*t1772*t2041*t2179;
  t3932 = t3901 + t3908;
  t3935 = -1.*t2030*t2765*t2041;
  t3957 = t2021*t3932;
  t3959 = t3935 + t3957;
  t3978 = -1.*t2021*t2030*t2041;
  t3982 = -1.*t2765*t3932;
  t4013 = t3978 + t3982;
  t4159 = t2478*t2063*t2041;
  t4179 = t1528*t1772*t2063*t2179;
  t4181 = t4159 + t4179;
  t4189 = t1528*t2030*t2765*t2063;
  t4197 = t2021*t4181;
  t4203 = t4189 + t4197;
  t4221 = t2021*t1528*t2030*t2063;
  t4226 = -1.*t2765*t4181;
  t4234 = t4221 + t4226;
  t4390 = -1.*t2478*t932*t2041;
  t4414 = -1.*t932*t1528*t1772*t2179;
  t4458 = t4390 + t4414;
  t4467 = -1.*t932*t1528*t2030*t2765;
  t4468 = t2021*t4458;
  t4469 = t4467 + t4468;
  t4480 = -1.*t2021*t932*t1528*t2030;
  t4481 = -1.*t2765*t4458;
  t4492 = t4480 + t4481;
  t4687 = -1.*t1528*t2030*t1772*t2765;
  t4696 = -1.*t2021*t1528*t2179;
  t4699 = t4687 + t4696;
  t4728 = t2021*t1528*t2030*t1772;
  t4734 = -1.*t1528*t2765*t2179;
  t4750 = t4728 + t4734;
  t4869 = -1.*t2030*t2063*t2041;
  t4885 = -1.*t932*t2179;
  t4887 = t4869 + t4885;
  t4934 = t1772*t2765*t4887;
  t4984 = t2021*t3481;
  t4997 = t4934 + t4984;
  t5002 = -1.*t2021*t1772*t4887;
  t5004 = t2765*t3481;
  t5034 = t5002 + t5004;
  t5134 = t2030*t2063;
  t5136 = t932*t2041*t2179;
  t5142 = t5134 + t5136;
  t5151 = t1772*t2765*t2237;
  t5157 = t2021*t5142;
  t5169 = t5151 + t5157;
  t5182 = -1.*t2021*t1772*t2237;
  t5184 = t2765*t5142;
  t5195 = t5182 + t5184;
  t5363 = -1.*t1772*t2041;
  t5369 = t2478*t1528*t2179;
  t5376 = t5363 + t5369;
  t5493 = t1528*t1772*t2063;
  t5496 = -1.*t2478*t3481;
  t5502 = t5493 + t5496;
  t2705 = -1.*t932*t1528*t1772;
  t5555 = -1.*t2478*t5142;
  t5569 = t2705 + t5555;
  t5411 = t2478*t2041;
  t5417 = t1528*t1772*t2179;
  t5418 = t5411 + t5417;
  t5723 = -1.*t1528*t2030*t2765;
  t5727 = -1.*t2021*t5418;
  t5731 = t5723 + t5727;
  t5761 = t2021*t1528*t2030;
  t5770 = -1.*t2765*t5418;
  t5805 = t5761 + t5770;
  t5900 = -1.*t2765*t3442;
  t5910 = -1.*t2021*t3520;
  t5925 = t5900 + t5910;
  t3824 = t2907*t3768;
  t5602 = t2478*t932*t1528;
  t5622 = -1.*t1772*t5142;
  t5629 = t5602 + t5622;
  t6032 = -1.*t932*t2030*t2041;
  t6034 = t2063*t2179;
  t6040 = t6032 + t6034;
  t6073 = -1.*t2765*t6040;
  t6079 = -1.*t2021*t5629;
  t6082 = t6073 + t6079;
  t6157 = t2021*t6040;
  t6162 = -1.*t2765*t5629;
  t6168 = t6157 + t6162;
  t6331 = t1528*t2030*t2765;
  t6333 = t2021*t5418;
  t6344 = t6331 + t6333;
  t5821 = t2907*t5805;
  t5870 = -1.*t2982*t5805;
  t6320 = 0.3455*t2907;
  t6327 = -0.35*t2982;
  t6330 = t6320 + t6327;
  t6374 = 0.35*t2907;
  t6376 = t6374 + t2992;
  t3814 = -1.*t2982*t3686;
  t3828 = t3814 + t3824;
  t6015 = -1.*t2982*t3768;
  t6537 = t2765*t6040;
  t6546 = t2021*t5629;
  t6548 = t6537 + t6546;
  t6252 = t2907*t6168;
  t6297 = -1.*t2982*t6168;
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
  p_output1[10]=0.3455*t2026*t2237 + 0.072*t2537*t2662 + 0.19875*(t2478*t2662 + t2705) + 0.3455*t2765*t2897 + t3036*t3110 + t3233*t3263 + 0.3455*(-1.*t2982*t3110 + t2907*t3263) - 0.7*(t2907*t3110 + t2982*t3263) + 0.072*t1528*t1772*t932;
  p_output1[11]=0.072*t1528*t1772*t2063 + 0.3455*t2026*t3442 + 0.072*t2537*t3481 + 0.19875*(-1.*t1528*t1772*t2063 + t2478*t3481) + 0.3455*t2765*t3520 + t3036*t3686 + t3233*t3768 - 0.7*(t2907*t3686 + t2982*t3768) + 0.3455*t3828;
  p_output1[12]=-0.072*t1528*t1772 - 0.3455*t2026*t2030*t2041 + 0.19875*(t1528*t1772 + t2041*t2179*t2478) + 0.072*t2041*t2179*t2537 + 0.3455*t2765*t3932 + t3036*t3959 + t3233*t4013 + 0.3455*(-1.*t2982*t3959 + t2907*t4013) - 0.7*(t2907*t3959 + t2982*t4013);
  p_output1[13]=0.3455*t1528*t2026*t2030*t2063 - 0.072*t1772*t2041*t2063 + 0.19875*(t1772*t2041*t2063 - 1.*t1528*t2063*t2179*t2478) - 0.072*t1528*t2063*t2179*t2537 + 0.3455*t2765*t4181 + t3036*t4203 + t3233*t4234 + 0.3455*(-1.*t2982*t4203 + t2907*t4234) - 0.7*(t2907*t4203 + t2982*t4234);
  p_output1[14]=0.3455*t2765*t4458 + t3036*t4469 + t3233*t4492 + 0.3455*(-1.*t2982*t4469 + t2907*t4492) - 0.7*(t2907*t4469 + t2982*t4492) - 0.3455*t1528*t2026*t2030*t932 + 0.072*t1772*t2041*t932 + 0.072*t1528*t2179*t2537*t932 + 0.19875*(-1.*t1772*t2041*t932 + t1528*t2179*t2478*t932);
  p_output1[15]=-0.3455*t1528*t2026*t2179 - 0.19875*t1528*t2030*t2478 - 0.072*t1528*t2030*t2537 + 0.3455*t1528*t1772*t2030*t2765 + t3233*t4699 + t3036*t4750 - 0.7*(t2982*t4699 + t2907*t4750) + 0.3455*(t2907*t4699 - 1.*t2982*t4750);
  p_output1[16]=0.3455*t2026*t3481 + 0.19875*t2478*t4887 + 0.072*t2537*t4887 - 0.3455*t1772*t2765*t4887 + t3233*t4997 + t3036*t5034 - 0.7*(t2982*t4997 + t2907*t5034) + 0.3455*(t2907*t4997 - 1.*t2982*t5034);
  p_output1[17]=0.19875*t2237*t2478 + 0.072*t2237*t2537 - 0.3455*t1772*t2237*t2765 + 0.3455*t2026*t5142 + t3233*t5169 + t3036*t5195 - 0.7*(t2982*t5169 + t2907*t5195) + 0.3455*(t2907*t5169 - 1.*t2982*t5195);
  p_output1[18]=0;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=0;
  p_output1[22]=0;
  p_output1[23]=0;
  p_output1[24]=0;
  p_output1[25]=0;
  p_output1[26]=0;
  p_output1[27]=-0.072*t1528*t1772*t2179 - 0.072*t2041*t2478 + 0.3455*t2765*t5376 + t2021*t3036*t5376 - 1.*t2765*t3233*t5376 + 0.3455*(-1.*t2765*t2907*t5376 - 1.*t2021*t2982*t5376) - 0.7*(t2021*t2907*t5376 - 1.*t2765*t2982*t5376) + 0.19875*t5418;
  p_output1[28]=0.072*t1528*t2063*t2478 + 0.072*t1772*t3481 + 0.19875*t3520 + 0.3455*t2765*t5502 + t2021*t3036*t5502 - 1.*t2765*t3233*t5502 + 0.3455*(-1.*t2765*t2907*t5502 - 1.*t2021*t2982*t5502) - 0.7*(t2021*t2907*t5502 - 1.*t2765*t2982*t5502);
  p_output1[29]=0.072*t1772*t5142 + 0.3455*t2765*t5569 + t2021*t3036*t5569 - 1.*t2765*t3233*t5569 + 0.3455*(-1.*t2765*t2907*t5569 - 1.*t2021*t2982*t5569) - 0.7*(t2021*t2907*t5569 - 1.*t2765*t2982*t5569) + 0.19875*t5629 - 0.072*t1528*t2478*t932;
  p_output1[30]=0.3455*t1528*t2030*t2765 + 0.3455*t2021*t5418 + t3233*t5731 + t3036*t5805 - 0.7*(t2982*t5731 + t5821) + 0.3455*(t2907*t5731 + t5870);
  p_output1[31]=0.3455*t2765*t3442 + 0.3455*t2021*t3520 + t3036*t3768 + t3233*t5925 - 0.7*(t3824 + t2982*t5925) + 0.3455*(t2907*t5925 + t6015);
  p_output1[32]=0.3455*t2021*t5629 + 0.3455*t2765*t6040 + t3233*t6082 + t3036*t6168 - 0.7*(t2982*t6082 + t6252) + 0.3455*(t2907*t6082 + t6297);
  p_output1[33]=t6330*t6344 + 0.3455*(t5870 - 1.*t2907*t6344) - 0.7*(t5821 - 1.*t2982*t6344) + t5805*t6376;
  p_output1[34]=-0.7*t3828 + 0.3455*(-1.*t2907*t3686 + t6015) + t3686*t6330 + t3768*t6376;
  p_output1[35]=t6168*t6376 + t6330*t6548 + 0.3455*(t6297 - 1.*t2907*t6548) - 0.7*(t6252 - 1.*t2982*t6548);
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

#include "J_FL.hh"

namespace SymFunction
{

void J_FL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
