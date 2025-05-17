#ifndef MHE_EST_HPP
#define MHE_EST_HPP

// #include <rclcpp/rclcpp.hpp>
#include <Eigen/Sparse>
#include <thread>

#include "decentral_legged_est/MheSrb.hpp"
#include "decentral_legged_est/EigenUtils.hpp"
#include "decentral_legged_est/Spline/Bezier_simple.hpp"

#include "../../src/Expressions/Mmat1_b1_description.hh"
#include "../../src/Expressions/Mmat2_b1_description.hh"
#include "../../src/Expressions/Mmat3_b1_description.hh"
#include "../../src/Expressions/Mmat4_b1_description.hh"
#include "../../src/Expressions/Mmat5_b1_description.hh"
#include "../../src/Expressions/Mmat6_b1_description.hh"
#include "../../src/Expressions/Mmat7_b1_description.hh"
#include "../../src/Expressions/Mmat8_b1_description.hh"
#include "../../src/Expressions/Mmat9_b1_description.hh"
#include "../../src/Expressions/Mmat10_b1_description.hh"
#include "../../src/Expressions/Mmat11_b1_description.hh"
#include "../../src/Expressions/Mmat12_b1_description.hh"
#include "../../src/Expressions/Mmat13_b1_description.hh"
#include "../../src/Expressions/Mmat14_b1_description.hh"
#include "../../src/Expressions/Ge_vec_b1_description.hh"
#include "../../src/Expressions/u_map_b1_description.hh"
#include "../../src/Expressions/C_mat_go1.hh"
#include "../../src/Expressions/J_FR.hh"
#include "../../src/Expressions/J_FL.hh"
#include "../../src/Expressions/J_RR.hh"
#include "../../src/Expressions/J_RL.hh"

using namespace Eigen;
// enum estimation_type {
//     MHE = 0,
//     KF = 1,
// };

struct robot_params
{
    // ros params
    std::vector<double> p_process_std_;
    std::vector<double> accel_input_std_;
    std::vector<double> accel_bias_std_;
    std::vector<double> gyro_input_std_;

    std::vector<double> GM_process_std_;
    std::vector<double> GRF_process_std_;
    std::vector<double> GM_Torso_process_std_;
    std::vector<double> Force_Torso_process_std_;

    std::vector<double> GM_meas_gyro_std_;
    std::vector<double> GM_meas_foot_std_;

    std::vector<double> quaternion_ib_;
    std::vector<double> p_ib_;

    int num_legs_;
    int num_joints_;
    int using_lo_v_;
    int using_lo_p_;
    int using_impact_map_;

    std::vector<double> joint_position_std_;
    std::vector<double> joint_velocity_std_;
    std::vector<double> foot_slide_std_;
    std::vector<double> foot_swing_std_;
    double contact_effort_theshold_;

    std::vector<double> p_init_std_;
    std::vector<double> v_init_std_;
    std::vector<double> foot_init_std_;
    std::vector<double> accel_bias_init_std_;
    std::vector<double> GM_init_std_;
    std::vector<double> GRF_init_std_;
    std::vector<double> GM_Torso_init_std_;
    std::vector<double> Force_Torso_init_std_;

    std::vector<double> vo_p_std_;

    // estimator params
    int rate_;
    int N_;
    int est_type_;
    int est_GRF_;
    int est_Torso_;

    // osqp params
    double rho_;
    double alpha_;
    double delta_;
    double sigma_;
    bool verbose_;
    bool adaptRho_;
    bool polish_;
    int maxQPIter_;
    double realtiveTol_;
    double absTol_;
    double primTol_;
    double dualTol_;
    double timeLimit_;
};

struct robot_store
{
    // IMU
    double imu_time_;
    Vector3d accel_b_;
    Vector3d angular_b_;

    // Encoder & Contact
    VectorXd joint_states_position_;
    VectorXd joint_states_velocity_;
    VectorXd joint_states_effort_;
    VectorXd contact_;

    VectorXd foot_anguler_rate_;

    MatrixXd p_imu_2_foot_;
    MatrixXd J_imu_2_foot_;

    // VO
    double vo_time_pre_;
    double vo_time_now_;
    bool vo_new_ = false;
    Vector3d vo_p_body_pre_2_body_;
    Quaterniond vo_quaternion_;

    // Decentralized Filter
    Quaterniond quaternion_;
    Quaterniond offset_quaternion_;

    Vector3d gt_p_;
    Vector3d gt_v_s_;
};

class DecentralizedEstimation
{
public:
    DecentralizedEstimation();

    void initialize(std::shared_ptr<robot_store> sub_ptr, std::shared_ptr<robot_params> params_ptr);
    void update(int T);
    void reset();
    VectorXd b_meas_;
    VectorXd contact_;

    VectorXd Disturbance;

    int solve_num = 1;
    double solve_time = 0.0;

private:
    double dt_;
    int N_;
    int est_type_;
    int est_GRF_;
    int est_Torso_;

    Vector3d gravity_;
    double contact_effort_theshold_ = 150.0;   // if using theshold to detect contact
    MatrixXd R_ib_ = Matrix3d::Identity(3, 3); // Go1 body frame is choosen to be Unitree_URDF_center frame during codegen
    MatrixXd p_ib_ = MatrixXd::Zero(1, 3);

private:
    std::shared_ptr<robot_store> robot_sub_ptr_; // robot instantaneous sensory data store
    std::shared_ptr<robot_params> params_ptr_;
    MHEproblem mhe_qp_;
    Bezier vo_curve_;

private:
    //---------------------------------------------------------------
    // Current Measurements
    Vector3d accel_b_;
    Vector3d accel_s_;
    Vector3d angular_b_;
    Matrix3d omega_skew_ = Matrix3d::Zero();

    Vector3d euler_sb_ = Vector3d::Zero();
    Vector3d euler_rate_sb_ = Vector3d::Zero();

    VectorXd joint_position_;
    VectorXd joint_velocity_;
    VectorXd joint_effort_;

    MatrixXd p_imu_2_foot_;
    MatrixXd J_imu_2_foot_;

    VectorXd foot_anguler_rate_;

    VectorXd contact_change_;

    Vector3d gt_v_s_;

    bool vo_new_meas_flag_ = false;
    double imu_time_;

    //---------------------------------------------------------------
    // Input && Measurement stack
    std::vector<int> discrete_time_stack;
    std::vector<double> imu_time_stack_;
    std::vector<Vector3d> accel_s_input_stack_;    // Spatial acceleration stacks; haven't correct for bias;
    std::vector<Vector3d> angular_b_input_stack_;  // Body gyro stacks; haven't correct for bias;
    std::vector<Matrix3d> R_input_rotation_stack_; // Spatial orientation stacks; R_sb, body frame in the world frame;
    std::vector<VectorXd> euler_input_rotation_stack_;

    std::vector<MatrixXd> p_imu_2_foot_stack_;
    std::vector<MatrixXd> J_imu_2_foot_stack_;
    std::vector<VectorXd> joint_position_stack_;
    std::vector<VectorXd> joint_velocity_stack_;
    std::vector<VectorXd> joint_effort_stack_;

    std::vector<VectorXd> contact_input_stack_;  // Contact stacks; 0, no contact, 1, contact
    std::vector<VectorXd> contact_change_stack_; // Contact stacks; 0, no contact, 1, contact

    std::vector<Vector3d> gt_v_s_input_stack_; // Body gyro stacks; haven't correct for bias;

    //---------------------------------------------------------------
    // VO params
    std::vector<int> vo_insert_idx_stack_;
    std::vector<int> vo_insert_discrete_time_stack_;

    Vector3d vo_p_body_pre_2_body_ = Vector3d::Zero();
    bool vo_to_be_processed_flag_ = false;
    Matrix3d R_vo_sb_pre_ = Matrix3d::Zero();

    // Debug stack
    // std::vector<Vector3d> vo_pose_body_pre_2_body_stack;
    // std::vector<Matrix3d> R_vo_sb_rotation_stack_; // Spatial orientation stacks; R_sb, body frame in the world frame when got camera;
    // std::vector<double> vo_time_stack_;
    // std::vector<int> vo_sychron_time_stack_;
private:
    // Contact Constraint
    int dim_contact_eq_;
    // SparseMatrix<double> A_contact_eq_;

    MatrixXd A_contact_eq_dense_;

    VectorXd b_contact_eq_;

    // Contact Constraint
    int dim_contact_ineq_;
    SparseMatrix<double> A_contact_ineq_;
    VectorXd b_contact_ineq_up_;
    VectorXd b_contact_ineq_lo_;

private:
    // Dynanmics
    MatrixXd Mmat_;
    MatrixXd Cmat_;
    MatrixXd Bmat_;
    MatrixXd Gvec_;
    VectorXd q_append_;
    VectorXd dq_append_;
    MatrixXd Cmat_T_p_;
    MatrixXd Cmat_T_Ralpha_;

    MatrixXd Mmat_p_;
    MatrixXd Mmat_p_Ralpha_;

    int using_impact_map_;

    MatrixXd Q_GM_meas_;
private:
    VectorXd GM_residual_;
    VectorXd GM_residual_integration_;
    VectorXd GM_;
    VectorXd GM_0_;
    MatrixXd K_;

private:
    int lo_p_start_idx_;
    int lo_v_start_idx_;
    int p_foot_start_idx_;
    int GM_torso_start_idx_;
    int wrench_torso_start_idx_;

private:
    //---------------------------------------------------------------
    // LO Measurement cost&constraints term
    // 0.5 * || v_i ||^2 _{Q_meas_}
    // A_meas x_i - b_meas - v_i = 0
    int dim_meas_;

    // int leg_odom_type_;
    int using_lo_v_;
    int using_lo_p_;

    int num_legs_;

    int dim_leg_;

    int dim_dof_;

    SparseMatrix<double> Identity_meas_;

    // A_meas_:
    // foot velocity: [  0   I   0;  ]
    // foot_position: [  -I   0   0  I;  ]
    SparseMatrix<double> A_meas_;

    // b_meas_:
    // foot_position: [  R_sb * fk   ];
    // foot_velocity: [  -R_sb * J * dq - R_sb omega.cross(fk)   ];
    // VectorXd b_meas_;

    // Q_meas_:
    //  foot_position: [   R_sb * (J * C_encoder_position_^{-1} * J')^{-1} * R_sb' ];
    //  foot_velocity: [   R_sb * ([-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]')^{-1} * R_sb' ];
    SparseMatrix<double> Q_meas_;
    MatrixXd C_meas_input_ = MatrixXd::Zero(9, 9);
    MatrixXd C_encoder_position_ = Matrix3d::Zero();
    MatrixXd C_encoder_velocity_ = Matrix3d::Zero();
    MatrixXd C_gyro_ = Matrix3d::Zero();
    MatrixXd C_foot_swing_ = Matrix3d::Zero();
    MatrixXd Q_foot_swing_ = Matrix3d::Zero();
    MatrixXd C_GM_omega_alpha_;

    //---------------------------------------------------------------
    // Camera cost&constraints term
    // 0.5 * || A_cam_pre_ * x_{T_pre} - A_cam_now_ * x_T - b_cam_ ||^2_{Q_cam_}
    // A_cam_pre * x_k - A_cam_now * x_k+1 - vcam_k = b_cam_k
    int dim_cam_;

    SparseMatrix<double> Identity_cam_meas_;

    // A_cam_pre:
    //   [  I  0   0   0;];
    SparseMatrix<double> A_cam_pre_;
    // A_cam_now:
    //   [  I  0   0   0;];
    SparseMatrix<double> A_cam_now_;

    // b_meas:
    // [    - R_sb_pre * vo_realtive_p  ;   ]
    VectorXd b_cam_;

    // Q_cam_pre:
    //  [   R_sb_pre * (Q_vo_p_)^{-1} * R_sb_pre'  ;  ]
    SparseMatrix<double> Q_cam_;
    MatrixXd Q_vo_p_ = Matrix3d::Zero();

private:
    //---------------------------------------------------------------
    // Dynamics cost terms
    // 0.5 * ||w_i||^2 _{Q_dyn}
    // A_dyn_ x_i - x_{i+1} - b_dyn -w_i = 0
    int dim_state_;

    SparseMatrix<double> Identity_dyn_;

    // A_dyn_:
    // [    I   dt* I   - 0.5 * dt^2 * R_sb     0;  ]
    // [    0   I       - dt * R_sb             0;  ]
    // [    0   0       I                       0;  ]
    // [    0   0       0                       I;  ] foot_position
    SparseMatrix<double> A_dyn_;

    // b_dyn_:
    // [    - 0.5 * dt^2 * accel_s; ]   p_s
    // [    - dt * accel_s        ; ]   v_s
    // [    0                     ; ]   accel_bias_b
    // [    0                     ; ]   foot_position
    VectorXd b_dyn_;

    // Q_dyn_:
    // (G_dyn * diag[C_p_, C_accel_, C_accel_bias_, C_foot_] * G_dyn').inverse()
    // G_dyn:
    // [  R_sb * dt          -0.5 * R_sb *dt^2   0;         0;  ]
    // [  0                  -R_sb * dt          0;         0;  ]
    // [  0                  0                   I * dt;    0;  ]
    // [  0                  0                   0;         R_sb * dt;  ] foot_position
    SparseMatrix<double> Q_dyn_;

    MatrixXd C_dyn_pv_ = MatrixXd::Zero(6, 6);

    MatrixXd C_p_ = Matrix3d::Zero();
    MatrixXd C_accel_ = Matrix3d::Zero();
    MatrixXd C_foot_slide_ = Matrix3d::Zero();
    MatrixXd Q_foot_slide_ = Matrix3d::Zero();
    MatrixXd C_accel_bias_ = Matrix3d::Zero();
    MatrixXd Q_accel_bias_ = Matrix3d::Zero();

    //---------------------------------------------------------------
    // GM & GRF dynamics
    MatrixXd C_GM_process_ = Matrix3d::Zero();
    MatrixXd C_GRF_process_ = Matrix3d::Zero();
    MatrixXd Q_GM_process_ = Matrix3d::Zero();
    MatrixXd Q_GRF_process_ = Matrix3d::Zero();

    MatrixXd Q_GM_Torso_process_ = MatrixXd::Zero(6, 6);
    MatrixXd Q_Force_Torso_process_ = MatrixXd::Zero(6, 6);

    MatrixXd C_GM_Torso_process_ = MatrixXd::Zero(6, 6);
    MatrixXd C_Force_Torso_process_ = MatrixXd::Zero(6, 6);

    // Block identity selection matrix
    MatrixXd S_;
    MatrixXd S_torso_;

    //---------------------------------------------------------------
    // Prior cost term
    // 0.5 * || x_0 - x_prior_0 ||^2 _{Q_prior_0}

    // x_prior:
    // [    0                       ;   ]   p_s
    // [    0                       ;   ]   v_s
    // [    0                       ;   ]   accel_bias_b
    // [    R_sb * fk               ;   ]   foot_position
    VectorXd x_prior_;

    // Q_prior_:
    //   [  Q_p_init  0          0                  0;   ]
    //   [  0         Q_v_init   0                  0;   ]
    //   [  0         0          Q_accel_bias_init  0;   ]
    //   [  0         0          0                  Q_foot_init;   ]  foot_position
    SparseMatrix<double> Q_prior_; // Gain prior

private:
    SparseMatrix<double> A_lcp_;

    VectorXd b_lcp_;

public:
    Matrix3d R_sb_;
    Vector3d p_vo_accmulate_ = Vector3d::Zero();

    //---------------------------------------------------------------
    // Moving horizon estimation
    VectorXd x_MHE_;
    Vector3d v_MHE_b_ = Vector3d::Zero();
    //---------------------------------------------------------------
    // Kalmanfilter
    VectorXd x_KF_;
    MatrixXd C_KF_;
    MatrixXd K_KF_;
    Vector3d v_KF_b_ = Vector3d::Zero();

private:
    void ConfigurateOSQP();

    void InitializeMHE();
    void UpdateMHE(int discrete_time);
    void GetMeasurement(int discrete_time);
    void UpdateVOConstraints(int discrete_time);

    void InitializeMomentumResidual();
    void updateMomentumResidual(int discrete_time);

    void InitializeKF();
    void UpdateKF();
    void MarginalizeKF();

    void Mmat_b1(MatrixXd &Mmat, const VectorXd &q);
    void StdVec2CovMat(const std::vector<double> &std, MatrixXd &Cov);   // Wrapper convert Vector3d std to Diag Covariance Matrix
    void StdVec2GainMat(const std::vector<double> &std, MatrixXd &Gain); // Wrapper convert Vector3d std to Diag Gain Matrix
    void tic(std::string str, int mode = 0);
    void toc(std::string str);
};
#endif // MHE_EST_HPP