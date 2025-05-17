/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:05:47 GMT-05:00
 */

#include "math2mat.hpp"
#include "mdefs.hpp"

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t3;
  double t4;
  double t8;
  double t5;
  double t15;
  double t11;
  double t18;
  double t38;
  double t39;
  double t40;
  double t37;
  double t41;
  double t42;
  double t44;
  double t45;
  double t46;
  double t50;
  double t51;
  double t52;
  double t6;
  double t12;
  double t60;
  double t67;
  double t16;
  double t17;
  double t19;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t25;
  double t26;
  double t28;
  double t29;
  double t30;
  double t32;
  double t33;
  double t34;
  double t43;
  double t47;
  double t48;
  double t49;
  double t53;
  double t54;
  double t55;
  double t56;
  double t57;
  double t58;
  double t61;
  double t62;
  double t63;
  double t64;
  double t65;
  double t66;
  double t68;
  double t69;
  double t70;
  double t72;
  double t73;
  double t74;
  double t75;
  double t76;
  double t77;
  double t78;
  double t79;
  double t80;
  double t81;
  double t86;
  double t87;
  double t88;
  double t89;
  double t90;
  double t91;
  double t92;
  double t27;
  double t31;
  double t35;
  double t36;
  double t105;
  double t106;
  double t107;
  double t108;
  double t59;
  double t71;
  double t82;
  double t83;
  double t109;
  double t110;
  double t111;
  double t112;
  double t127;
  double t128;
  double t129;
  double t130;
  double t143;
  double t144;
  double t145;
  double t84;
  double t85;
  double t93;
  double t94;
  double t113;
  double t114;
  double t115;
  double t116;
  double t131;
  double t132;
  double t133;
  double t134;
  double t138;
  double t139;
  double t140;
  double t141;
  double t148;
  double t149;
  double t150;
  double t151;
  double t172;
  double t173;
  double t174;
  double t168;
  double t169;
  double t170;
  double t164;
  double t165;
  double t166;
  double t95;
  double t96;
  double t97;
  double t117;
  double t118;
  double t119;
  double t135;
  double t136;
  double t137;
  double t176;
  double t177;
  double t178;
  double t189;
  double t190;
  double t191;
  t3 = Cos(var1[4]);
  t4 = Power(t3,2);
  t8 = Sin(var1[4]);
  t5 = Cos(var1[5]);
  t15 = Sin(var1[3]);
  t11 = Sin(var1[5]);
  t18 = Cos(var1[3]);
  t38 = 0.00875*t5;
  t39 = -0.002184*t11;
  t40 = t38 + t39;
  t37 = 0.002934*t8;
  t41 = t3*t40;
  t42 = t37 + t41;
  t44 = 0.002934*t3;
  t45 = -1.*t8*t40;
  t46 = t44 + t45;
  t50 = 0.002184*t5;
  t51 = 0.00875*t11;
  t52 = t50 + t51;
  t6 = Power(t5,2);
  t12 = Power(t11,2);
  t60 = t5*t52;
  t67 = -1.*t52*t11;
  t16 = -30.247*t3*t15*t8;
  t17 = t5*t15*t8;
  t19 = t18*t11;
  t20 = t17 + t19;
  t21 = 30.247*t3*t5*t20;
  t22 = t18*t5;
  t23 = -1.*t15*t8*t11;
  t24 = t22 + t23;
  t25 = -30.247*t3*t11*t24;
  t26 = t16 + t21 + t25;
  t28 = -1.*t18*t5*t8;
  t29 = t15*t11;
  t30 = t28 + t29;
  t32 = t5*t15;
  t33 = t18*t8*t11;
  t34 = t32 + t33;
  t43 = t8*t42;
  t47 = t3*t46;
  t48 = t43 + t47;
  t49 = -1.*t3*t48*t11;
  t53 = -1.*t5*t52;
  t54 = t3*t42*t11;
  t55 = -1.*t8*t46*t11;
  t56 = t53 + t54 + t55;
  t57 = t8*t56;
  t58 = t49 + t57;
  t61 = -1.*t3*t42*t11;
  t62 = t8*t46*t11;
  t63 = t60 + t61 + t62;
  t64 = t3*t5*t63;
  t65 = -1.*t3*t5*t42;
  t66 = t5*t8*t46;
  t68 = t65 + t66 + t67;
  t69 = -1.*t3*t11*t68;
  t70 = t64 + t69;
  t72 = -1.*t8*t42;
  t73 = -1.*t3*t46;
  t74 = t72 + t73;
  t75 = t3*t5*t74;
  t76 = t3*t5*t42;
  t77 = -1.*t5*t8*t46;
  t78 = t52*t11;
  t79 = t76 + t77 + t78;
  t80 = t8*t79;
  t81 = t75 + t80;
  t86 = -1.*t40*t11;
  t87 = t60 + t86;
  t88 = t11*t87;
  t89 = -1.*t5*t40;
  t90 = t89 + t67;
  t91 = t5*t90;
  t92 = t88 + t91;
  t27 = 30.247*t18*t3*t8;
  t31 = 30.247*t3*t5*t30;
  t35 = -30.247*t3*t11*t34;
  t36 = t27 + t31 + t35;
  t105 = -30.247*t18*t4*t15;
  t106 = 30.247*t20*t30;
  t107 = 30.247*t34*t24;
  t108 = t105 + t106 + t107;
  t59 = 30.247*t3*t5*t58;
  t71 = 30.247*t8*t70;
  t82 = -30.247*t3*t11*t81;
  t83 = t59 + t71 + t82;
  t109 = 30.247*t20*t58;
  t110 = -30.247*t3*t15*t70;
  t111 = 30.247*t24*t81;
  t112 = t109 + t110 + t111;
  t127 = 30.247*t30*t58;
  t128 = 30.247*t18*t3*t70;
  t129 = 30.247*t34*t81;
  t130 = t127 + t128 + t129;
  t143 = -0.027977*t3*t5;
  t144 = 0.81937*t8;
  t145 = -0.000189*t3*t11;
  t84 = 0.088744698*t3*t6;
  t85 = 0.088744698*t3*t12;
  t93 = 30.247*t8*t92;
  t94 = t84 + t85 + t93;
  t113 = 0.088744698*t5*t20;
  t114 = -0.088744698*t11*t24;
  t115 = -30.247*t3*t15*t92;
  t116 = t113 + t114 + t115;
  t131 = 0.088744698*t5*t30;
  t132 = -0.088744698*t11*t34;
  t133 = 30.247*t18*t3*t92;
  t134 = t131 + t132 + t133;
  t138 = -0.001395*t3*t5;
  t139 = 0.000189*t8;
  t140 = -0.787698*t3*t11;
  t141 = t138 + t139 + t140;
  t148 = 0.188947*t3*t5;
  t149 = -0.027977*t8;
  t150 = 0.001395*t3*t11;
  t151 = t148 + t149 + t150;
  t172 = 0.088744698*t5*t58;
  t173 = 30.247*t70*t92;
  t174 = -0.088744698*t11*t81;
  t168 = 0.787698*t5;
  t169 = -0.001395*t11;
  t170 = t168 + t169;
  t164 = -0.001395*t5;
  t165 = 0.188947*t11;
  t166 = t164 + t165;
  t95 = -0.066059448*t3*t5;
  t96 = -0.26466125*t3*t11;
  t97 = t95 + t96;
  t117 = -0.066059448*t20;
  t118 = 0.26466125*t24;
  t119 = t117 + t118;
  t135 = -0.066059448*t30;
  t136 = 0.26466125*t34;
  t137 = t135 + t136;
  t176 = -0.066059448*t58;
  t177 = 0.26466125*t81;
  t178 = t143 + t144 + t145 + t176 + t177;
  t189 = -4.818420431999999e-6*t5;
  t190 = -0.0287535161075*t11;
  t191 = t189 + t190;
  p_output1[0]=30.247*t12*t4 + 30.247*t4*t6 + 30.247*Power(t8,2);
  p_output1[1]=t26;
  p_output1[2]=t36;
  p_output1[3]=t83;
  p_output1[4]=t94;
  p_output1[5]=t97;
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=0;
  p_output1[15]=0;
  p_output1[16]=0;
  p_output1[17]=0;
  p_output1[18]=t26;
  p_output1[19]=30.247*Power(t20,2) + 30.247*Power(t24,2) + 30.247*Power(t15,2)*t4;
  p_output1[20]=t108;
  p_output1[21]=t112;
  p_output1[22]=t116;
  p_output1[23]=t119;
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
  p_output1[36]=t36;
  p_output1[37]=t108;
  p_output1[38]=30.247*Power(t30,2) + 30.247*Power(t34,2) + 30.247*Power(t18,2)*t4;
  p_output1[39]=t130;
  p_output1[40]=t134;
  p_output1[41]=t137;
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
  p_output1[54]=t83;
  p_output1[55]=t112;
  p_output1[56]=t130;
  p_output1[57]=-1.*t11*t141*t3 + t151*t3*t5 + 30.247*Power(t58,2) + 30.247*Power(t70,2) + (t143 + t144 + t145)*t8 + 30.247*Power(t81,2);
  p_output1[58]=t172 + t173 + t174 - 1.*t11*t170*t3 + t166*t3*t5 + (-0.027977*t11 + 0.000189*t5)*t8;
  p_output1[59]=t178;
  p_output1[60]=0;
  p_output1[61]=0;
  p_output1[62]=0;
  p_output1[63]=0;
  p_output1[64]=0;
  p_output1[65]=0;
  p_output1[66]=0;
  p_output1[67]=0;
  p_output1[68]=0;
  p_output1[69]=0;
  p_output1[70]=0;
  p_output1[71]=0;
  p_output1[72]=t94;
  p_output1[73]=t116;
  p_output1[74]=t134;
  p_output1[75]=t11*t151 + t172 + t173 + t174 + t141*t5;
  p_output1[76]=0.000260376943932*t12 + t11*t166 + t170*t5 + 0.000260376943932*t6 + 30.247*Power(t92,2);
  p_output1[77]=t191;
  p_output1[78]=0;
  p_output1[79]=0;
  p_output1[80]=0;
  p_output1[81]=0;
  p_output1[82]=0;
  p_output1[83]=0;
  p_output1[84]=0;
  p_output1[85]=0;
  p_output1[86]=0;
  p_output1[87]=0;
  p_output1[88]=0;
  p_output1[89]=0;
  p_output1[90]=t97;
  p_output1[91]=t119;
  p_output1[92]=t137;
  p_output1[93]=t178;
  p_output1[94]=t191;
  p_output1[95]=0.821830059771932;
  p_output1[96]=0;
  p_output1[97]=0;
  p_output1[98]=0;
  p_output1[99]=0;
  p_output1[100]=0;
  p_output1[101]=0;
  p_output1[102]=0;
  p_output1[103]=0;
  p_output1[104]=0;
  p_output1[105]=0;
  p_output1[106]=0;
  p_output1[107]=0;
  p_output1[108]=0;
  p_output1[109]=0;
  p_output1[110]=0;
  p_output1[111]=0;
  p_output1[112]=0;
  p_output1[113]=0;
  p_output1[114]=0;
  p_output1[115]=0;
  p_output1[116]=0;
  p_output1[117]=0;
  p_output1[118]=0;
  p_output1[119]=0;
  p_output1[120]=0;
  p_output1[121]=0;
  p_output1[122]=0;
  p_output1[123]=0;
  p_output1[124]=0;
  p_output1[125]=0;
  p_output1[126]=0;
  p_output1[127]=0;
  p_output1[128]=0;
  p_output1[129]=0;
  p_output1[130]=0;
  p_output1[131]=0;
  p_output1[132]=0;
  p_output1[133]=0;
  p_output1[134]=0;
  p_output1[135]=0;
  p_output1[136]=0;
  p_output1[137]=0;
  p_output1[138]=0;
  p_output1[139]=0;
  p_output1[140]=0;
  p_output1[141]=0;
  p_output1[142]=0;
  p_output1[143]=0;
  p_output1[144]=0;
  p_output1[145]=0;
  p_output1[146]=0;
  p_output1[147]=0;
  p_output1[148]=0;
  p_output1[149]=0;
  p_output1[150]=0;
  p_output1[151]=0;
  p_output1[152]=0;
  p_output1[153]=0;
  p_output1[154]=0;
  p_output1[155]=0;
  p_output1[156]=0;
  p_output1[157]=0;
  p_output1[158]=0;
  p_output1[159]=0;
  p_output1[160]=0;
  p_output1[161]=0;
  p_output1[162]=0;
  p_output1[163]=0;
  p_output1[164]=0;
  p_output1[165]=0;
  p_output1[166]=0;
  p_output1[167]=0;
  p_output1[168]=0;
  p_output1[169]=0;
  p_output1[170]=0;
  p_output1[171]=0;
  p_output1[172]=0;
  p_output1[173]=0;
  p_output1[174]=0;
  p_output1[175]=0;
  p_output1[176]=0;
  p_output1[177]=0;
  p_output1[178]=0;
  p_output1[179]=0;
  p_output1[180]=0;
  p_output1[181]=0;
  p_output1[182]=0;
  p_output1[183]=0;
  p_output1[184]=0;
  p_output1[185]=0;
  p_output1[186]=0;
  p_output1[187]=0;
  p_output1[188]=0;
  p_output1[189]=0;
  p_output1[190]=0;
  p_output1[191]=0;
  p_output1[192]=0;
  p_output1[193]=0;
  p_output1[194]=0;
  p_output1[195]=0;
  p_output1[196]=0;
  p_output1[197]=0;
  p_output1[198]=0;
  p_output1[199]=0;
  p_output1[200]=0;
  p_output1[201]=0;
  p_output1[202]=0;
  p_output1[203]=0;
  p_output1[204]=0;
  p_output1[205]=0;
  p_output1[206]=0;
  p_output1[207]=0;
  p_output1[208]=0;
  p_output1[209]=0;
  p_output1[210]=0;
  p_output1[211]=0;
  p_output1[212]=0;
  p_output1[213]=0;
  p_output1[214]=0;
  p_output1[215]=0;
  p_output1[216]=0;
  p_output1[217]=0;
  p_output1[218]=0;
  p_output1[219]=0;
  p_output1[220]=0;
  p_output1[221]=0;
  p_output1[222]=0;
  p_output1[223]=0;
  p_output1[224]=0;
  p_output1[225]=0;
  p_output1[226]=0;
  p_output1[227]=0;
  p_output1[228]=0;
  p_output1[229]=0;
  p_output1[230]=0;
  p_output1[231]=0;
  p_output1[232]=0;
  p_output1[233]=0;
  p_output1[234]=0;
  p_output1[235]=0;
  p_output1[236]=0;
  p_output1[237]=0;
  p_output1[238]=0;
  p_output1[239]=0;
  p_output1[240]=0;
  p_output1[241]=0;
  p_output1[242]=0;
  p_output1[243]=0;
  p_output1[244]=0;
  p_output1[245]=0;
  p_output1[246]=0;
  p_output1[247]=0;
  p_output1[248]=0;
  p_output1[249]=0;
  p_output1[250]=0;
  p_output1[251]=0;
  p_output1[252]=0;
  p_output1[253]=0;
  p_output1[254]=0;
  p_output1[255]=0;
  p_output1[256]=0;
  p_output1[257]=0;
  p_output1[258]=0;
  p_output1[259]=0;
  p_output1[260]=0;
  p_output1[261]=0;
  p_output1[262]=0;
  p_output1[263]=0;
  p_output1[264]=0;
  p_output1[265]=0;
  p_output1[266]=0;
  p_output1[267]=0;
  p_output1[268]=0;
  p_output1[269]=0;
  p_output1[270]=0;
  p_output1[271]=0;
  p_output1[272]=0;
  p_output1[273]=0;
  p_output1[274]=0;
  p_output1[275]=0;
  p_output1[276]=0;
  p_output1[277]=0;
  p_output1[278]=0;
  p_output1[279]=0;
  p_output1[280]=0;
  p_output1[281]=0;
  p_output1[282]=0;
  p_output1[283]=0;
  p_output1[284]=0;
  p_output1[285]=0;
  p_output1[286]=0;
  p_output1[287]=0;
  p_output1[288]=0;
  p_output1[289]=0;
  p_output1[290]=0;
  p_output1[291]=0;
  p_output1[292]=0;
  p_output1[293]=0;
  p_output1[294]=0;
  p_output1[295]=0;
  p_output1[296]=0;
  p_output1[297]=0;
  p_output1[298]=0;
  p_output1[299]=0;
  p_output1[300]=0;
  p_output1[301]=0;
  p_output1[302]=0;
  p_output1[303]=0;
  p_output1[304]=0;
  p_output1[305]=0;
  p_output1[306]=0;
  p_output1[307]=0;
  p_output1[308]=0;
  p_output1[309]=0;
  p_output1[310]=0;
  p_output1[311]=0;
  p_output1[312]=0;
  p_output1[313]=0;
  p_output1[314]=0;
  p_output1[315]=0;
  p_output1[316]=0;
  p_output1[317]=0;
  p_output1[318]=0;
  p_output1[319]=0;
  p_output1[320]=0;
  p_output1[321]=0;
  p_output1[322]=0;
  p_output1[323]=0;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 18, (mwSize) 18, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "Mmat1_b1_description.hh"

namespace SymFunction
{

void Mmat1_b1_description_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
