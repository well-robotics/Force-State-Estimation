/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:43 GMT-05:00
 */

#ifndef C_MAT_GO1_HH
#define C_MAT_GO1_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void C_mat_go1_raw(double *p_output1, const double *var1,const double *var2);

  inline void C_mat_go1(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 18, 1);
    assert_size_matrix(var2, 18, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 18, 18);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    C_mat_go1_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // C_MAT_GO1_HH
