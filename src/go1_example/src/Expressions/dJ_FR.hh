/*
 * Automatically Generated from Mathematica.
 * Sun 22 Sep 2024 09:12:33 GMT-05:00
 */

#ifndef DJ_FR_HH
#define DJ_FR_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void dJ_FR_raw(double *p_output1, const double *var1);

  inline void dJ_FR(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 36, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 3, 18);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    dJ_FR_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // DJ_FR_HH
