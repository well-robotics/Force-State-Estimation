/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:28 GMT-05:00
 */

#ifndef U_VEC_B1_DESCRIPTION_HH
#define U_VEC_B1_DESCRIPTION_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void u_vec_b1_description_raw(double *p_output1, const double *var1,const double *var2);

  inline void u_vec_b1_description(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 18, 1);
    assert_size_matrix(var2, 12, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 18, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    u_vec_b1_description_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // U_VEC_B1_DESCRIPTION_HH
