/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:05:54 GMT-05:00
 */

#ifndef MMAT11_B1_DESCRIPTION_HH
#define MMAT11_B1_DESCRIPTION_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Mmat11_b1_description_raw(double *p_output1, const double *var1);

  inline void Mmat11_b1_description(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 18, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 18, 18);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Mmat11_b1_description_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // MMAT11_B1_DESCRIPTION_HH
