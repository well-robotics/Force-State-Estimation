#include "decentral_legged_est/DecentralEst.hpp"
#include <fstream>

// Function pointer type
typedef void (*MmatFunction)(MatrixXd &, const VectorXd &);

// Array of function pointers
MmatFunction MmatFunctions[] = {
    SymFunction::Mmat1_b1_description,
    SymFunction::Mmat2_b1_description,
    SymFunction::Mmat3_b1_description,
    SymFunction::Mmat4_b1_description,
    SymFunction::Mmat5_b1_description,
    SymFunction::Mmat6_b1_description,
    SymFunction::Mmat7_b1_description,
    SymFunction::Mmat8_b1_description,
    SymFunction::Mmat9_b1_description,
    SymFunction::Mmat10_b1_description,
    SymFunction::Mmat11_b1_description,
    SymFunction::Mmat12_b1_description,
    SymFunction::Mmat13_b1_description,
    SymFunction::Mmat14_b1_description,
};

// Function pointer type
typedef void (*JFunction)(MatrixXd &, const VectorXd &);

// Array of function pointers
JFunction JFunctions[] = {
    SymFunction::J_FR,
    SymFunction::J_FL,
    SymFunction::J_RR,
    SymFunction::J_RL,
};

const int numFunctions = sizeof(MmatFunctions) / sizeof(MmatFunctions[0]);

DecentralizedEstimation::DecentralizedEstimation()
{
}

// Set the Estimation parameters and QP solver
void DecentralizedEstimation::initialize(std::shared_ptr<robot_store> sub, std::shared_ptr<robot_params> params)
{
    robot_sub_ptr_ = sub;
    params_ptr_ = params;

    est_type_ = params_ptr_->est_type_;
    est_GRF_ = params_ptr_->est_GRF_;
    est_Torso_ = params_ptr_->est_Torso_;

    dt_ = 1.0 / params_ptr_->rate_;
    N_ = params_ptr_->N_;

    num_legs_ = params_ptr_->num_legs_;

    using_lo_v_ = params_ptr_->using_lo_v_;
    using_lo_p_ = params_ptr_->using_lo_p_;
    using_impact_map_ = params_ptr_->using_impact_map_;

    dim_leg_ = params_ptr_->num_joints_;
    dim_dof_ = 6 + dim_leg_ * num_legs_;

    dim_state_ = 9 + est_GRF_ * (6 + dim_leg_ * num_legs_ + num_legs_ * 3) + est_GRF_ * est_Torso_ * 6 + using_lo_p_ * 3 * num_legs_;
    dim_meas_ = est_GRF_ * (6 + dim_leg_ * num_legs_) + (using_lo_v_ + using_lo_p_) * 3 * num_legs_;

    dim_cam_ = 3;
    dim_contact_eq_ = num_legs_ * 3;
    dim_contact_ineq_ = est_GRF_ * 1 * num_legs_;

    mhe_qp_.setHorizon(N_, dim_state_, dim_meas_, dim_cam_);

    contact_effort_theshold_ = params_ptr_->contact_effort_theshold_;
    gravity_ << 0, 0, -9.81;

    Quaterniond q_ib;
    q_ib.w() = params_ptr_->quaternion_ib_[0];
    q_ib.x() = params_ptr_->quaternion_ib_[1];
    q_ib.y() = params_ptr_->quaternion_ib_[2];
    q_ib.z() = params_ptr_->quaternion_ib_[3];
    R_ib_ = q_ib.toRotationMatrix();
    p_ib_ << params_ptr_->p_ib_[0], params_ptr_->p_ib_[1], params_ptr_->p_ib_[2];

    //---------------------------------------------------------------
    // Gain setup from params; Q: gains, inv of covariance; C: covariance
    StdVec2CovMat(params_ptr_->p_process_std_, C_p_);
    StdVec2CovMat(params_ptr_->accel_input_std_, C_accel_);
    StdVec2CovMat(params_ptr_->accel_bias_std_, C_accel_bias_);
    StdVec2CovMat(params_ptr_->joint_position_std_, C_encoder_position_);
    StdVec2CovMat(params_ptr_->joint_velocity_std_, C_encoder_velocity_);
    StdVec2CovMat(params_ptr_->gyro_input_std_, C_gyro_);
    StdVec2CovMat(params_ptr_->foot_slide_std_, C_foot_slide_);
    StdVec2CovMat(params_ptr_->foot_swing_std_, C_foot_swing_);

    StdVec2GainMat(params_ptr_->accel_bias_std_, Q_accel_bias_);
    StdVec2GainMat(params_ptr_->foot_slide_std_, Q_foot_slide_);
    StdVec2GainMat(params_ptr_->foot_swing_std_, Q_foot_swing_);
    StdVec2GainMat(params_ptr_->vo_p_std_, Q_vo_p_);

    std::vector<double> GM_std;
    GM_std.push_back(params_ptr_->GM_meas_gyro_std_[0]);
    GM_std.push_back(params_ptr_->GM_meas_gyro_std_[1]);
    GM_std.push_back(params_ptr_->GM_meas_gyro_std_[2]);
    GM_std.push_back(params_ptr_->GM_meas_gyro_std_[3]);
    GM_std.push_back(params_ptr_->GM_meas_gyro_std_[4]);
    GM_std.push_back(params_ptr_->GM_meas_gyro_std_[5]);

    for (int i = 0; i < num_legs_; i++)
    {
        for (int j = 0; j < dim_leg_; j++)
        {
            GM_std.push_back(params_ptr_->GM_meas_foot_std_[j]);
        }
    }
    Q_GM_meas_ = MatrixXd::Zero(dim_dof_, dim_dof_);
    StdVec2GainMat(GM_std, Q_GM_meas_);
    StdVec2CovMat(params_ptr_->GM_process_std_, C_GM_process_);
    StdVec2CovMat(params_ptr_->GRF_process_std_, C_GRF_process_);
    StdVec2GainMat(params_ptr_->GM_process_std_, Q_GM_process_);
    StdVec2GainMat(params_ptr_->GRF_process_std_, Q_GRF_process_);

    StdVec2GainMat(params_ptr_->GM_Torso_process_std_, Q_GM_Torso_process_);
    StdVec2GainMat(params_ptr_->Force_Torso_process_std_, Q_Force_Torso_process_);

    StdVec2CovMat(params_ptr_->GM_Torso_process_std_, C_GM_Torso_process_);
    StdVec2CovMat(params_ptr_->Force_Torso_process_std_, C_Force_Torso_process_);

    // Initialized the matrix size
    x_MHE_ = VectorXd::Zero(dim_state_);
    x_KF_ = VectorXd::Zero(dim_state_);

    // Declare the sparse matrix to the proper size
    Q_prior_.resize(dim_state_, dim_state_);
    Q_prior_.setZero();
    x_prior_ = VectorXd::Zero(dim_state_);

    Identity_dyn_.resize(dim_state_, dim_state_);
    Identity_dyn_.setIdentity();
    A_dyn_.resize(dim_state_, dim_state_);
    A_dyn_.setIdentity();
    Q_dyn_.resize(dim_state_, dim_state_);
    Q_dyn_.setZero();
    b_dyn_ = VectorXd::Zero(dim_state_);
    C_dyn_pv_.block<3, 3>(0, 0) = C_p_;
    C_dyn_pv_.block<3, 3>(3, 3) = C_accel_;

    Identity_meas_.resize(dim_meas_, dim_meas_);
    Identity_meas_.setIdentity();
    A_meas_.resize(dim_meas_, dim_state_);
    A_meas_.setZero();
    Q_meas_.resize(dim_meas_, dim_meas_);
    Q_meas_.setZero();
    b_meas_ = VectorXd::Zero(dim_meas_);
    C_meas_input_.block<3, 3>(0, 0) = C_encoder_velocity_;
    C_meas_input_.block<3, 3>(3, 3) = C_encoder_position_;
    C_meas_input_.block<3, 3>(6, 6) = C_gyro_;

    Identity_cam_meas_.resize(dim_cam_, dim_cam_);
    Identity_cam_meas_.setIdentity();
    A_cam_pre_.resize(dim_cam_, dim_state_);
    A_cam_pre_.setZero();
    A_cam_now_.resize(dim_cam_, dim_state_);
    A_cam_now_.setZero();
    Q_cam_.resize(dim_cam_, dim_cam_);
    Q_cam_.setZero();
    b_cam_ = VectorXd::Zero(dim_cam_);

    // A_contact_eq_.resize(dim_contact_eq_, dim_state_);
    // A_contact_eq_.setZero();

    A_contact_eq_dense_ = MatrixXd::Zero(dim_contact_eq_, dim_state_);
    b_contact_eq_ = VectorXd::Zero(dim_contact_eq_);

    A_contact_ineq_.resize(dim_contact_ineq_, dim_state_);
    A_contact_ineq_.setZero();
    b_contact_ineq_up_ = VectorXd::Zero(dim_contact_ineq_);
    b_contact_ineq_lo_ = VectorXd::Zero(dim_contact_ineq_);

    // Dynamics
    q_append_ = VectorXd::Zero(dim_dof_);
    dq_append_ = VectorXd::Zero(dim_dof_);

    // Mmat_ = [Mmat_p_  Mmat_p_Ralpha_]
    Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
    Mmat_p_ = MatrixXd::Zero(dim_dof_, 3);
    Mmat_p_Ralpha_ = MatrixXd::Zero(dim_dof_, dim_dof_ - 3);

    // Cmat_T_ = [Cmat_T_p_  Cmat_T_Ralpha_]
    Cmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
    Cmat_T_p_ = MatrixXd::Zero(dim_dof_, 3);
    Cmat_T_Ralpha_ = MatrixXd::Zero(dim_dof_, dim_dof_ - 3);

    Bmat_ = MatrixXd::Zero(dim_dof_, dim_leg_ * num_legs_);
    Gvec_ = MatrixXd::Zero(dim_dof_, 1);

    //---------------------------------------------------------------
    // Const matrix setup
    // Selection matrix setup
    S_ = MatrixXd::Zero(dim_leg_ * num_legs_, dim_dof_);

    for (int i = 0; i < num_legs_; i++)
    {
        S_.block(i * dim_leg_, 6 + i * dim_leg_, dim_leg_, dim_leg_) = MatrixXd::Identity(dim_leg_, dim_leg_);
    }

    // Selection Torso Matrix setup
    S_torso_ = MatrixXd::Zero(6, dim_dof_);
    S_torso_.block(0, 0, 6, 6) = MatrixXd::Identity(6, 6);

    // A_cam_pre_/A_cam_now_ * x:
    //  [  I   0   0;  ]
    A_cam_pre_.coeffRef(0, 0) = 1.0;
    A_cam_pre_.coeffRef(1, 1) = 1.0;
    A_cam_pre_.coeffRef(2, 2) = 1.0;
    A_cam_now_.coeffRef(0, 0) = 1.0;
    A_cam_now_.coeffRef(1, 1) = 1.0;
    A_cam_now_.coeffRef(2, 2) = 1.0;

    //---------------------------------------------------------------
    // Config indexing
    lo_v_start_idx_ = est_GRF_ * (num_legs_ * dim_leg_ + 6);
    lo_p_start_idx_ = est_GRF_ * (num_legs_ * dim_leg_ + 6) +
                      using_lo_v_ * 3 * num_legs_;

    p_foot_start_idx_ = 9 +
                        est_GRF_ * (6 + num_legs_ * dim_leg_ + num_legs_ * 3) +
                        est_Torso_ * 6;

    GM_torso_start_idx_ = 9;
    wrench_torso_start_idx_ = 9 + 6 + num_legs_ * dim_leg_ + 3 * num_legs_;

    Disturbance = VectorXd::Zero(dim_dof_);

    switch (est_type_)
    {
    case 0:
    {
        InitializeMHE();
        break;
    }
    case 1:
    {
        InitializeKF();
        UpdateKF();
        break;
    }
    case 2:
    {
        InitializeMomentumResidual();
    }
    default:
    {
        std::cout << est_type_ + " not a valid estimation type." << std::endl;
        break;
    }
    }
}

void DecentralizedEstimation::update(int T)
{
    switch (est_type_)
    {
    case 0:
    {
        UpdateMHE(T);
        // std::cout << "updateMHE to " + std::to_string(T) << std::endl;

        if (vo_to_be_processed_flag_)
        {
            UpdateVOConstraints(T);
            vo_to_be_processed_flag_ = false;
        }

        if (T >= N_)
        {
            mhe_qp_.marginalizeQP(T - N_); // the marginalization is independent of update() and assemble_image()
        }
        tic("formulate and solve");
        mhe_qp_.initQP(T); // update the osqp, note to keep the sparsity structure
        // toc("formulate and solve");

        // tic("formulate and solve");
        mhe_qp_.solveQP();
        toc("formulate and solve");

        mhe_qp_.getsolution(T);

        x_MHE_ = mhe_qp_.x_solution_now;

        Vector3d p_imu_2_opti;
        // p_imu_2_opti << 0.0, 0.0, 0.0;
        // p_imu_2_opti << 0.016041, 0.089061, 0.0579875;
        v_MHE_b_ = R_sb_.transpose() * (x_MHE_.segment<3>(3) + angular_b_.cross(p_imu_2_opti));

        break;
    }
    case 1:
    {
        UpdateKF();
        Vector3d p_imu_2_opti;
        // p_imu_2_opti << 0.0, 0.0, 0.0;
        // p_imu_2_opti << 0.016041, 0.089061, 0.0579875;
        v_KF_b_ = R_sb_.transpose() * (x_KF_.segment<3>(3) + angular_b_.cross(p_imu_2_opti));
        break;
        std::cout << x_KF_ << std::endl;
    }
    case 2:
    {
        updateMomentumResidual(T);
        break;
    }
    default:
    {
        std::cout << est_type_ + " not a valid estimation type." << std::endl;
        break;
    }
    }
}

void DecentralizedEstimation::InitializeMomentumResidual()
{
    GM_residual_ = VectorXd::Zero(dim_dof_);
    GM_residual_integration_ = VectorXd::Zero(dim_dof_);
    Disturbance = VectorXd::Zero(dim_dof_);
    K_ = 80 * MatrixXd::Identity(dim_dof_, dim_dof_);
    GetMeasurement(0);

    Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
    Mmat_b1(Mmat_, q_append_); // To be edit to configurate euler

    Cmat_.setZero();

    dq_append_.segment<3>(0) = gt_v_s_;

    SymFunction::C_mat_go1(Cmat_, q_append_, dq_append_);
    MatrixXd Cmat_T = Cmat_.transpose();

    Bmat_ = MatrixXd::Zero(dim_dof_, dim_leg_ * num_legs_);
    SymFunction::u_map_b1_description(Bmat_, q_append_);

    Gvec_ = MatrixXd::Zero(dim_dof_, 1);
    SymFunction::Ge_vec_b1_description(Gvec_, q_append_, VectorXd::Zero(dim_dof_));
    Gvec_ = -Gvec_;

    GM_0_ = Mmat_ * dq_append_;

    GM_residual_integration_ += (Cmat_T * dq_append_ - Gvec_ + Bmat_ * joint_effort_ + GM_residual_) * dt_;

    GM_residual_ = K_ * (-GM_residual_integration_);

    if (est_Torso_)
    {
        Disturbance.segment<6>(0) = GM_residual_.segment<6>(0);

        for (int i = 0; i < num_legs_; i++)
        {
            MatrixXd S_i = S_.block(i * dim_leg_, 0, dim_leg_, dim_dof_);
            MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
            JFunctions[i](J_i, q_append_);
            MatrixXd J_i_T = J_i.transpose();
            Disturbance.segment<3>(6 + 3 * i) = (S_i * J_i_T).inverse() * S_i * GM_residual_;
            // if (contact_input_stack_.back()(i) == 1.0)
            // {
            //     Disturbance.segment<6>(0) = Disturbance.segment<6>(0) - S_torso_ * J_i_T * Disturbance.segment<3>(6 + 3 * i);
            // }
            Disturbance.segment<6>(0) = Disturbance.segment<6>(0) - S_torso_ * J_i_T * Disturbance.segment<3>(6 + 3 * i);
        }
    }
    else
    {

        MatrixXd J_append = MatrixXd::Zero(num_legs_ * 3, dim_dof_);

        for (int i = 0; i < num_legs_; i++)
        {
            MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
            JFunctions[i](J_i, q_append_);
            J_append.block(3 * i, 0, 3, dim_dof_) = J_i;
        }

        Disturbance = (J_append.transpose()).inverse() * GM_residual_;
    }
}

void DecentralizedEstimation::updateMomentumResidual(int T)
{
    Disturbance = VectorXd::Zero(dim_dof_);

    // int impact_count = (contact_change_stack_.back().array() == 1.0).count();
    // MatrixXd J_impact_B = MatrixXd::Zero(dim_dof_, 3 * impact_count);
    // MatrixXd J_impact_C = MatrixXd::Zero(3 * impact_count, dim_dof_);
    // int foot_cout = 0;
    // MatrixXd Mmat_inv = Mmat_.inverse();

    // for (int i = 0; i < num_legs_; i++)
    // {
    //     if (contact_change_stack_.back()(i) == 1.0)
    //     {
    //         MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
    //         JFunctions[i](J_i, q_append_);
    //         J_impact_B.block(0, 0 + foot_cout * 3, dim_dof_, 3) = -J_i.transpose();
    //         J_impact_C.block(0 + foot_cout * 3, 0, 3, dim_dof_) = -J_i * Mmat_inv;
    //         foot_cout += 1;
    //     }
    // }

    // MatrixXd A_impact = MatrixXd::Identity(dim_dof_, dim_dof_) - J_impact_B * (J_impact_C * J_impact_B).inverse() * J_impact_C;

    GetMeasurement(T);

    Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
    Mmat_b1(Mmat_, q_append_); // To be edit to configurate euler

    Cmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);

    dq_append_.segment<3>(0) = gt_v_s_;
    // dq_append_.segment<3>(3) = angular_b_in    Cmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);

    SymFunction::Ge_vec_b1_description(Gvec_, q_append_, VectorXd::Zero(dim_dof_));
    Gvec_ = -Gvec_;
    GM_ = Mmat_ * dq_append_;

    GM_residual_integration_ += (Cmat_ * dq_append_ - Gvec_ + Bmat_ * joint_effort_ + GM_residual_) * dt_;

    GM_residual_ = K_ * (GM_ - GM_residual_integration_ - GM_0_);
    Disturbance.segment<6>(0) = GM_residual_.segment<6>(0);
    for (int i = 0; i < num_legs_; i++)
    {
        MatrixXd S_i = S_.block(i * dim_leg_, 0, dim_leg_, dim_dof_);
        MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
        JFunctions[i](J_i, q_append_);
        MatrixXd J_i_T = J_i.transpose();
        Disturbance.segment<3>(6 + 3 * i) = (S_i * J_i_T).inverse() * S_i * GM_residual_;
        // if (contact_input_stack_.back()(i) == 1.0)
        // {
        //     Disturbance.segment<6>(0) = Disturbance.segment<6>(0) - S_torso_ * J_i_T * Disturbance.segment<3>(6 + 3 * i);
        // }
        Disturbance.segment<6>(0) = Disturbance.segment<6>(0) - S_torso_ * J_i_T * Disturbance.segment<3>(6 + 3 * i);
    }
    std::cout << Disturbance << std::endl;
}

void DecentralizedEstimation::ConfigurateOSQP()
{
    //---------------------------------------------------------------
    // OSQP configurate
    mhe_qp_.osqp.settings()->setWarmStart(false);
    mhe_qp_.osqp.settings()->setAdaptiveRho(params_ptr_->adaptRho_);
    mhe_qp_.osqp.settings()->setVerbosity(params_ptr_->verbose_);
    mhe_qp_.osqp.settings()->setPolish(params_ptr_->polish_);
    mhe_qp_.osqp.settings()->setMaxIteration(params_ptr_->maxQPIter_);
    mhe_qp_.osqp.settings()->setRho(params_ptr_->rho_);
    mhe_qp_.osqp.settings()->setAlpha(params_ptr_->alpha_);
    mhe_qp_.osqp.settings()->setDelta(params_ptr_->delta_);
    mhe_qp_.osqp.settings()->setSigma(params_ptr_->sigma_);
    mhe_qp_.osqp.settings()->setRelativeTolerance(params_ptr_->realtiveTol_);
    mhe_qp_.osqp.settings()->setAbsoluteTolerance(params_ptr_->absTol_);
    mhe_qp_.osqp.settings()->setPrimalInfeasibilityTolerance(params_ptr_->primTol_);
    mhe_qp_.osqp.settings()->setDualInfeasibilityTolerance(params_ptr_->dualTol_);
    mhe_qp_.osqp.settings()->setTimeLimit(params_ptr_->timeLimit_);
}

void DecentralizedEstimation::InitializeMHE()
{

    ConfigurateOSQP();

    GetMeasurement(0);

    MatrixXd Q_p_int = Matrix3d::Zero();
    MatrixXd Q_v_int = Matrix3d::Zero();
    MatrixXd Q_accel_bias_init = Matrix3d::Zero();
    MatrixXd Q_GM_init = Matrix3d::Zero();
    MatrixXd Q_GRF_init = Matrix3d::Zero();
    MatrixXd Q_GM_Torso_init = MatrixXd::Zero(6, 6);
    MatrixXd Q_Force_Torso_init = MatrixXd::Zero(6, 6);

    StdVec2GainMat(params_ptr_->p_init_std_, Q_p_int);
    StdVec2GainMat(params_ptr_->v_init_std_, Q_v_int);
    StdVec2GainMat(params_ptr_->accel_bias_init_std_, Q_accel_bias_init);
    StdVec2GainMat(params_ptr_->GM_init_std_, Q_GM_init);
    StdVec2GainMat(params_ptr_->GRF_init_std_, Q_GRF_init);
    StdVec2GainMat(params_ptr_->GM_Torso_init_std_, Q_GM_Torso_init);
    StdVec2GainMat(params_ptr_->Force_Torso_init_std_, Q_Force_Torso_init);

    //---------------------------------------------------------------
    // Prior cost at 0
    // 0.5 * || x_0 - x_prior ||^2_{Q_prior_}

    // x_prior:
    // [    0                       ;   ]   p_s
    // [    0                       ;   ]   v_s
    // [    0                       ;   ]   accel_bias_b
    // [    R_sb * fk               ;   ]   foot_position
    x_prior_.segment<3>(0) << 0.0, 0.0, 0.0;
    x_prior_.segment<3>(3) << 0.0, 0.0, 0.0;
    x_prior_.segment<3>(6) << 0.0, 0.0, 0.0;

    // Q_prior_:
    //   [  Q_p_init  0          0                  0;   ]
    //   [  0         Q_v_init   0                  0;   ]
    //   [  0         0          Q_accel_bias_init  0;   ]
    //   [  0         0          0                  Q_foot_init;   ]  foot_position
    Q_prior_.setZero();
    Q_meas_.setZero();

    EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_, 0, 0, Q_p_int);
    EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_, 3, 3, Q_v_int);
    EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_, 6, 6, Q_accel_bias_init);

    //---------------------------------------------------------------
    // Meas cost & constraints at 0
    // 0.5 * || A_meas_ * x_0 - b_meas_ ||^2_{Q_meas_}
    A_meas_.setZero();
    MatrixXd Q_meas_dense = MatrixXd::Zero(dim_meas_, dim_meas_);

    //---------------------------------------------------------------
    // GM measurement at 0
    if (est_GRF_)
    {
        Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
        Mmat_b1(Mmat_, q_append_);

        Mmat_p_ = Mmat_.block(0, 0, dim_dof_, 3);
        Mmat_p_Ralpha_ = Mmat_.block(0, 3, dim_dof_, dim_dof_ - 3);
        // std::cout << Mmat_ << std::endl;
        //      p   v              accel_b_   GM_i   GRF_i
        //   [  0   - Mmat_p_      0          I      0;  ]
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    3,
                                                    -Mmat_p_); // v
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    9,
                                                    MatrixXd::Identity(dim_dof_, dim_dof_)); // GM
        b_meas_.segment(0, dim_dof_) = Mmat_p_Ralpha_ * dq_append_.segment(3, dim_dof_ - 3);
        // Q_meas_dense.block(0, 0, dim_dof_, dim_dof_) = (Mmat_p_Ralpha_ *
        //                                                 C_GM_omega_alpha_ *
        //                                                 Mmat_p_Ralpha_.transpose())
        //                                                    .inverse();

        Q_meas_dense.block(0, 0, dim_dof_, dim_dof_) = Q_GM_meas_;
        // prior
        x_prior_.segment(9, dim_dof_) << b_meas_.segment(0, dim_dof_);
        EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_,
                                                    9,
                                                    9,
                                                    Q_GM_Torso_init);
        for (int i = 0; i < num_legs_; i++)
        {
            x_prior_.segment(9 + 6 + num_legs_ * dim_leg_ + i * 3, 3) << 0.0, 0.0, 0.0; // GRF_foot_prior
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_,
                                                        9 + 6 + i * dim_leg_,
                                                        9 + 6 + i * dim_leg_,
                                                        Q_GM_init);
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_,
                                                        9 + 6 + num_legs_ * dim_leg_ + i * 3,
                                                        9 + 6 + num_legs_ * dim_leg_ + i * 3,
                                                        Q_GRF_init);
        }

        if (est_Torso_)
        {
            // prior
            x_prior_.segment(9 + 6 + num_legs_ * (dim_leg_ + 3), 6) << VectorXd::Zero(6); // Force_Torso_prior
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_,
                                                        9 + 6 + num_legs_ * (dim_leg_ + 3),
                                                        9 + 6 + num_legs_ * (dim_leg_ + 3),
                                                        Q_Force_Torso_init);
        }
    }
    //---------------------------------------------------------------
    // LO measurement at 0
    // Q_meas_:
    //  foot_position: [   R_sb * (J * C_encoder_position_^{-1} * J')^{-1} * R_sb' ];
    //  foot_velocity: [   R_sb * ([-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]')^{-1} * R_sb' ];

    // b_meas_:
    //  foot_position: [  R_sb * fk   ];
    //  foot_velocity: [  -R_sb * J * dq - R_sb omega.cross(fk)   ];
    switch (using_lo_v_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {
        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_v_start_idx_ + i * 3,
                                                        3,
                                                        MatrixXd::Identity(3, 3)); // v

            b_meas_.segment<3>(lo_v_start_idx_ + i * 3) = -R_sb_ *
                                                              J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                                                              joint_velocity_.segment<3>(i * 3) -
                                                          R_sb_ *
                                                              angular_b_.cross(
                                                                  p_imu_2_foot_.block<3, 1>(i * 3, 0));
            if (contact_(i) == 0.0)
            {
                Q_meas_dense.block<3, 3>(lo_v_start_idx_ + i * 3, lo_v_start_idx_ + i * 3) << Q_foot_swing_;
            }
            else
            {
                MatrixXd G_meas_i = MatrixXd::Zero(3, 9);

                MatrixXd C_meas_dense = MatrixXd::Zero(3, 3);

                G_meas_i.block<3, 3>(0, 0) = -J_imu_2_foot_.block<3, 3>(i * 3, 0);

                G_meas_i.block<3, 3>(0, 3) = -omega_skew_ * J_imu_2_foot_.block<3, 3>(i * 3, 0);

                Matrix3d kin_skew = Matrix3d::Zero();
                EigenUtils::vector3dSkew(kin_skew, p_imu_2_foot_.block<3, 1>(i * 3, 0));
                G_meas_i.block<3, 3>(0, 6) = kin_skew;

                C_meas_dense = R_sb_ * G_meas_i * C_meas_input_ *
                               G_meas_i.transpose() * R_sb_.transpose();

                Q_meas_dense.block<3, 3>(lo_v_start_idx_ + i * 3, lo_v_start_idx_ + i * 3) = C_meas_dense.inverse();
            }
        }
        break;
    }
    }

    switch (using_lo_p_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {

        MatrixXd Q_foot_init = Matrix3d::Zero();
        StdVec2GainMat(params_ptr_->foot_init_std_, Q_foot_init);

        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        0,
                                                        -MatrixXd::Identity(3, 3)); // p
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        p_foot_start_idx_ + i * 3,
                                                        MatrixXd::Identity(3, 3)); // p_foot

            b_meas_.segment<3>(lo_p_start_idx_ + i * 3) = R_sb_ *
                                                          p_imu_2_foot_.block<3, 1>(i * 3, 0);

            Q_meas_dense.block<3, 3>(lo_p_start_idx_ + i * 3,
                                     lo_p_start_idx_ + i * 3)
                << R_sb_ *
                       (J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                        C_encoder_position_ *
                        J_imu_2_foot_.block<3, 3>(i * 3, 0).transpose())
                           .inverse() *
                       R_sb_.transpose();

            // prior
            x_prior_.segment(p_foot_start_idx_ + i * 3, 3) = b_meas_.segment(lo_p_start_idx_ + i * 3, 3);
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_prior_,
                                                        p_foot_start_idx_ + i * 3,
                                                        p_foot_start_idx_ + i * 3,
                                                        Q_foot_init);
        }

        break;
    }
    default:
    {
        std::cout << est_type_ + " not a valid leg odom type." << std::endl;
        break;
    }
    }

    EigenUtils::SparseMatrixBlockAsignFromDense(Q_meas_, 0, 0, Q_meas_dense);

    // 0.5 * || x_0 - x_prior ||^2_{Q_prior_}
    mhe_qp_.addVariable("x_0", dim_state_);
    mhe_qp_.addCost("Prior_0", x_prior_, Q_prior_);
    mhe_qp_.addCostDependency("Prior_0", "x_0", Identity_dyn_);

    // // A_meas_ * x_i - v_{i} = b_meas_
    // mhe_qp_.addVariable("v_0", dim_meas_);
    // mhe_qp_.addConstraints("Measurement_0", b_meas_, b_meas_);
    // mhe_qp_.addConstraintDependency("Measurement_0", "x_0", A_meas_);
    // mhe_qp_.addConstraintDependency("Measurement_0", "v_0", -Identity_meas_);

    // // || v_i ||^2_{Q_meas_}
    // mhe_qp_.addCost("Measurement_0", VectorXd::Zero(dim_meas_), Q_meas_);
    // mhe_qp_.addCostDependency("Measurement_0", "v_0", Identity_meas_);

    // std::cout << Q_meas_dense << std::endl;
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(Q_meas_dense, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // std::cout << "rank: " << svd.rank() << std::endl;

    // || A_meas_ * x_i - b_meas_ ||^2_{Q_meas_}
    mhe_qp_.addCost("Measurement_0", b_meas_, Q_meas_);
    mhe_qp_.addCostDependency("Measurement_0", "x_0", A_meas_);

    //---------------------------------------------------------------
    // Contact Constraints
    if (est_GRF_ && 1)
    {
        // dim_contact_eq_ = static_cast<int>(num_legs_ - contact_.sum()) * 3;
        // dim_contact_eq_ = num_legs_ * 3;
        // A_contact_eq_.resize(dim_contact_eq_, dim_state_);
        // A_contact_eq_.setZero();

        A_contact_eq_dense_.setZero();
        b_contact_eq_ = VectorXd::Zero(dim_contact_eq_);
        // int id = 0;
        // A_contact_ineq_.setZero();
        for (int i = 0; i < num_legs_; i++)
        {
            //---------------------------------------------------------------
            // Contact Equality Constraints
            if (contact_(i) == 0.0)
            {
                // EigenUtils::SparseMatrixBlockAsignFromDense(A_contact_eq_, i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3, MatrixXd::Identity(3, 3)); // GRF
                A_contact_eq_dense_.block(i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3, 3, 3) = MatrixXd::Identity(3, 3); // GRF
                b_contact_eq_.segment<3>(i * 3) = VectorXd::Zero(3);
                // id++;
                // b_contact_ineq_lo_.segment<1>(i * 1) = -mhe_qp_.osqp_infinity * VectorXd::Ones(1);
                // b_contact_ineq_up_.segment<1>(i * 1) = mhe_qp_.osqp_infinity * VectorXd::Ones(1);
            }
            else
            {
                // EigenUtils::SparseMatrixBlockAsignFromDense(A_contact_eq_, i * 3, 9 + num_legs_ * dim_leg_ + i * 3, MatrixXd::Zero(3, 3)); // GRF

                A_contact_eq_dense_.block(i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3, 3, 3) = MatrixXd::Zero(3, 3); // GRF
                b_contact_eq_.segment<3>(i * 3) = VectorXd::Zero(3);
                // b_contact_ineq_lo_.segment<1>(i * 1) = -mhe_qp_.osqp_infinity * VectorXd::Ones(1);
                // b_contact_ineq_up_.segment<1>(i * 1) = mhe_qp_.osqp_infinity * VectorXd::Ones(1);
            }

            //---------------------------------------------------------------
            // Contact Inequality Constraints
            // EigenUtils::SparseMatrixBlockAsignFromDense(A_contact_ineq_, i * 1, 9 + num_legs_ * dim_leg_ + i * 3 + 2, MatrixXd::Identity(1, 1)); // GRF
        }

        mhe_qp_.addConstraints("Contact_0", b_contact_eq_, b_contact_eq_);
        mhe_qp_.addConstraintDependency("Contact_0", "x_0", A_contact_eq_dense_.sparseView());

        // mhe_qp_.addConstraints("Lcp_0", b_contact_ineq_lo_, b_contact_ineq_up_);
        // mhe_qp_.addConstraintDependency("Lcp_0", "x_0", A_contact_ineq_);
    }

    mhe_qp_.updateQP(0);
}

void DecentralizedEstimation::UpdateMHE(int T)
{
    std::string T_pre_string = std::to_string(T - 1);
    std::string T_string = std::to_string(T);
    std::string X_T_pre_string = "x_" + T_pre_string;
    std::string X_T_string = "x_" + T_string;

    std::string W_T_pre_string = "w_" + T_pre_string;
    std::string Vcam_T_pre_string = "vcam_" + T_pre_string;
    std::string V_T_string = "v_" + T_string;

    std::string Dyn_string = "Dynamic_" + T_pre_string;
    std::string Cam_Meas_string = "VO_measurement_" + T_pre_string;
    std::string Meas_string = "Measurement_" + T_string;
    std::string Contact_string = "Contact_" + T_string;
    // std::string Lcp_string = "Lcp_" + T_string;

    mhe_qp_.addVariable(W_T_pre_string, dim_state_);
    mhe_qp_.addVariable(Vcam_T_pre_string, dim_cam_);
    mhe_qp_.addVariable(X_T_string, dim_state_);
    // mhe_qp_.addVariable(V_T_string, dim_meas_);

    double dt = dt_;

    //---------------------------------------------------------------
    // Dynamic cost & constraints at T-1
    // 0.5 * || x_T - f( x_{T-1},u_{T-1} ) ||^2_{Q_{T-1}}
    // => 0.5 * || A_dyn_ x_{T-1} - x_T - b_dyn_ ||^2_{Q_dyn_{T-1}}

    // b_dyn_:
    // [    - 0.5 * dt_^2 * accel_s_; ]   p_s
    // [    - dt_ * accel_s_        ; ]   v_s
    // [    0                     ; ]   accel_bias_b
    // [    0                     ; ]   foot_position
    b_dyn_.segment(0, 3) << -dt_ * dt_ / 2 * accel_s_;
    b_dyn_.segment(3, 3) << -dt_ * accel_s_;

    // A_dyn_:
    // [    I   dt_* I   - 0.5 * dt_^2 * R_sb     0;  ]
    // [    0   I       - dt_ * R_sb             0;  ]
    // [    0   0       I                       0;  ]
    // [    0   0       0                       I;  ] foot_position
    A_dyn_.setIdentity();
    EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_, 0, 3, dt_ * Matrix3d::Identity());
    EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_, 0, 6, -dt_ * dt_ / 2 * R_sb_);
    EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_, 3, 6, -dt_ * R_sb_);

    // Q_dyn_:
    // (G_dyn * diag[C_p_, C_accel_, C_accel_bias_, C_foot_] * G_dyn').inverse()
    // G_dyn:
    // [  R_sb * dt_          -0.5 * R_sb *dt_^2   0;         0;  ]
    // [  0                  -R_sb * dt_          0;         0;  ]
    // [  0                  0                   I * dt_;    0;  ]
    // [  0                  0                   0;         R_sb * dt_;  ] foot_position
    Q_dyn_.setZero();

    MatrixXd G_dyn_pv = MatrixXd::Zero(6, 6);
    G_dyn_pv.block<3, 3>(0, 0) = R_sb_ * dt_;
    G_dyn_pv.block<3, 3>(0, 3) = 0.5 * R_sb_ * dt_ * dt_;
    G_dyn_pv.block<3, 3>(3, 3) = R_sb_ * dt_;
    MatrixXd Q_dyn_pv = (G_dyn_pv * C_dyn_pv_ * G_dyn_pv.transpose()).inverse();

    EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_, 0, 0, Q_dyn_pv);
    EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                6, 6,
                                                1 / (dt_ * dt_) * Q_accel_bias_);

    //---------------------------------------------------------------
    // GM, GRF states update
    // GM_i^+ = GM_i + dt_* S_i_ * C^T_p * v + dt_ * S_i_ * J_i * GRF_i + dt_ * {S_i_ * C^T_{w,alpha} - S_i_ * G  + B * u}
    // GRF_i^+ = GRF_i
    if (est_GRF_)
    {

        dq_append_.segment<3>(0) = x_MHE_.segment<3>(3);

        SymFunction::C_mat_go1(Cmat_, q_append_, dq_append_);
        Cmat_T_p_ = (Cmat_.block(0, 0, 3, dim_dof_)).transpose();
        Cmat_T_Ralpha_ = (Cmat_.block(3, 0, dim_dof_ - 3, dim_dof_)).transpose();

        SymFunction::u_map_b1_description(Bmat_, q_append_);

        SymFunction::Ge_vec_b1_description(Gvec_, q_append_, VectorXd::Zero(dim_dof_));
        Gvec_ = -Gvec_;

        int GM_start_idx = 9;
        int GRF_foot_start_idx = 9 + 6 + num_legs_ * dim_leg_;

        // A_dyn_:
        //      p   v                 accel_b_     GM_i     GRF_i
        //   [  0   dt *  C^T_p       0            I        dt  * J_i^T;  ]
        EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_,
                                                    GM_start_idx,
                                                    3,
                                                    dt * Cmat_T_p_); // v
        // b_dyn_:
        // - dt_ * {C^T_{pitch,alpha} [w;dalpha] - G  + B * u}
        b_dyn_.segment(GM_start_idx, dim_dof_) = -dt_ *
                                                 (Cmat_T_Ralpha_ *
                                                      dq_append_.segment(3, dim_dof_ - 3) -
                                                  Gvec_ + Bmat_ * joint_effort_);
        EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                    GM_start_idx,
                                                    GM_start_idx,
                                                    1 / (dt_ * dt_) * Q_GM_Torso_process_);

        for (int i = 0; i < num_legs_; i++)
        {

            MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
            JFunctions[i](J_i, q_append_);
            EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_,
                                                        GM_start_idx,
                                                        GRF_foot_start_idx + i * 3,
                                                        dt * J_i.transpose());
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                        GM_start_idx + 6 + i * dim_leg_,
                                                        GM_start_idx + 6 + i * dim_leg_,
                                                        1 / (dt_ * dt_) * Q_GM_process_);
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                        GRF_foot_start_idx + i * 3,
                                                        GRF_foot_start_idx + i * 3,
                                                        1 / (dt_ * dt_) * Q_GRF_process_);
        }

        if (est_Torso_)
        {

            MatrixXd J_torso = MatrixXd::Zero(6, dim_dof_);
            J_torso.block(0, 0, 6, 6) = MatrixXd::Identity(6, 6);
            EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_,
                                                        GM_start_idx, wrench_torso_start_idx_,
                                                        dt * J_torso.transpose()); // Wrench
            EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                        wrench_torso_start_idx_,
                                                        wrench_torso_start_idx_,
                                                        1 / (dt_ * dt_) * Q_Force_Torso_process_);
        }
    }

    //---------------------------------------------------------------
    // Leg states update

    switch (using_lo_p_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {

        for (int i = 0; i < num_legs_; i++)
        {
            if (contact_input_stack_.back()(i)) // if contact at foot_idx, fixed the stance foot position using large Q_foot_slide_
            {
                EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                            p_foot_start_idx_ + i * 3,
                                                            p_foot_start_idx_ + i * 3,
                                                            1 / (dt_ * dt_) *
                                                                R_sb_ * Q_foot_slide_ * R_sb_.transpose());
            }
            else // if no contact at foot_idx, relax the swing foot position using large Q_foot_swing_
            {
                EigenUtils::SparseMatrixBlockAsignFromDense(Q_dyn_,
                                                            p_foot_start_idx_ + i * 3,
                                                            p_foot_start_idx_ + i * 3,
                                                            1 / (dt_ * dt_) *
                                                                R_sb_ * Q_foot_swing_ * R_sb_.transpose());
            }
        }
        break;
    }
    default:
    {
        std::cout << using_lo_p_ + " not a valid leg odom type." << std::endl;
        break;
    }
    }

    // // A_dyn_ x_i - x_{i+1} - w_{i} = b_dyn_
    // mhe_qp_.addConstraints(Dyn_string, b_dyn_, b_dyn_);
    // mhe_qp_.addConstraintDependency(Dyn_string, W_T_pre_string, -Identity_dyn_);
    // mhe_qp_.addConstraintDependency(Dyn_string, X_T_string, -Identity_dyn_);
    // mhe_qp_.addConstraintDependency(Dyn_string, X_T_pre_string, A_dyn_);

    // // 0.5 * || w_{i} ||^2 _{Q_dyn_}
    // mhe_qp_.addCost(Dyn_string, VectorXd::Zero(dim_state_), Q_dyn_);
    // mhe_qp_.addCostDependency(Dyn_string, W_T_pre_string, Identity_dyn_);

    // Note: if no interpolation, this part could be better formulated
    //---------------------------------------------------------------
    // Camera_measurement cost & constraints at T
    // 0.5 * || A_cam_pre_ * x_{T_pre} - A_cam_now_ * x_T - b_cam_ ||^2_{Q_cam_}
    b_cam_.segment(0, 3) << mhe_qp_.osqp_infinity * Vector3d::Ones(); // unconstrained place holder, reserve for camera ahead of vo arrival

    Q_cam_.setZero();
    Matrix3d Q_cam_dense = R_sb_ * Q_vo_p_ * R_sb_.transpose();
    EigenUtils::SparseMatrixBlockAsignFromDense(Q_cam_, 0, 0, Q_cam_dense);

    // // A_cam_pre_ * x_pre - A_cam_now_ * x_now - vcam_pre = b_cam_pre; for consequetive states
    // mhe_qp_.addConstraints(Cam_Meas_string, -b_cam_, b_cam_);
    // mhe_qp_.addConstraintDependency(Cam_Meas_string, X_T_pre_string, A_cam_pre_);
    // mhe_qp_.addConstraintDependency(Cam_Meas_string, X_T_string, -A_cam_now_);
    // mhe_qp_.addConstraintDependency(Cam_Meas_string, Vcam_T_pre_string, -Identity_cam_meas_);

    // // 0.5 * || vcam_pre ||^2 _{Q_cam_pre}
    // mhe_qp_.addCost(Cam_Meas_string, VectorXd::Zero(dim_cam_), Q_cam_);
    // mhe_qp_.addCostDependency(Cam_Meas_string, Vcam_T_pre_string, Identity_cam_meas_);

    GetMeasurement(T);

    //---------------------------------------------------------------
    // Impact correction
    if (est_GRF_ && using_impact_map_)
    {
        int impact_count = (contact_change_stack_.back().array() == 1.0).count();
        MatrixXd J_impact_B = MatrixXd::Zero(dim_dof_, 3 * impact_count);
        MatrixXd J_impact_C = MatrixXd::Zero(3 * impact_count, dim_dof_);
        int foot_cout = 0;
        MatrixXd Mmat_inv = Mmat_.inverse();

        for (int i = 0; i < num_legs_; i++)
        {
            if (contact_change_stack_.back()(i) == 1.0)
            {
                MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
                JFunctions[i](J_i, q_append_);
                J_impact_B.block(0, 0 + foot_cout * 3, dim_dof_, 3) = -J_i.transpose();
                J_impact_C.block(0 + foot_cout * 3, 0, 3, dim_dof_) = J_i * Mmat_inv;
                foot_cout += 1;
            }
        }

        MatrixXd A_impact = MatrixXd::Identity(dim_dof_, dim_dof_) - J_impact_B * (J_impact_C * J_impact_B).inverse() * J_impact_C;
        std::cout << A_impact << std::endl;
        SparseMatrix<double> A_impact_reorder;
        A_impact_reorder.resize(dim_state_, dim_state_);
        A_impact_reorder.setIdentity();
        EigenUtils::SparseMatrixBlockReplaceFromDense(A_impact_reorder,
                                                      9,
                                                      9,
                                                      A_impact);
        A_dyn_ = A_dyn_ * A_impact_reorder;

        // MatrixXd F_map = (J_impact_C * J_impact_B).inverse() * J_impact_C;
        // VectorXd dq_local = dq_append_;
        // dq_local.segment<3>(0) = x_KF_.segment<3>(3);
        // std::cout << F_map * Mmat_ * dq_local << std::endl;
    }

    // A_dyn_ x_i - x_{i+1} - w_{i} = b_dyn_
    mhe_qp_.addConstraints(Dyn_string, b_dyn_, b_dyn_);
    mhe_qp_.addConstraintDependency(Dyn_string, W_T_pre_string, -Identity_dyn_);
    mhe_qp_.addConstraintDependency(Dyn_string, X_T_string, -Identity_dyn_);
    mhe_qp_.addConstraintDependency(Dyn_string, X_T_pre_string, A_dyn_);

    // 0.5 * || w_{i} ||^2 _{Q_dyn_}
    mhe_qp_.addCost(Dyn_string, VectorXd::Zero(dim_state_), Q_dyn_);
    mhe_qp_.addCostDependency(Dyn_string, W_T_pre_string, Identity_dyn_);

    // A_cam_pre_ * x_pre - A_cam_now_ * x_now - vcam_pre = b_cam_pre; for consequetive states
    mhe_qp_.addConstraints(Cam_Meas_string, -b_cam_, b_cam_);
    mhe_qp_.addConstraintDependency(Cam_Meas_string, X_T_pre_string, A_cam_pre_);
    mhe_qp_.addConstraintDependency(Cam_Meas_string, X_T_string, -A_cam_now_);
    mhe_qp_.addConstraintDependency(Cam_Meas_string, Vcam_T_pre_string, -Identity_cam_meas_);

    // 0.5 * || vcam_pre ||^2 _{Q_cam_pre}
    mhe_qp_.addCost(Cam_Meas_string, VectorXd::Zero(dim_cam_), Q_cam_);
    mhe_qp_.addCostDependency(Cam_Meas_string, Vcam_T_pre_string, Identity_cam_meas_);

    //---------------------------------------------------------------
    // Meas cost & constraints at T
    // 0.5 * || A_meas_ * x_T - b_meas_ ||^2_{Q_meas_}

    Q_meas_.setZero();
    MatrixXd Q_meas_dense = MatrixXd::Zero(dim_meas_, dim_meas_);

    A_meas_.setZero();

    //---------------------------------------------------------------
    // GM measurement at T
    if (est_GRF_)
    {
        Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
        Mmat_b1(Mmat_, q_append_);

        Mmat_p_ = Mmat_.block(0, 0, dim_dof_, 3);
        Mmat_p_Ralpha_ = Mmat_.block(0, 3, dim_dof_, dim_dof_ - 3);

        //      p   v              accel_b_   GM_i   GRF_i
        //   [  0   - Mmat_p_      0          I      0;  ]
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    3,
                                                    -Mmat_p_); // v
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    9,
                                                    MatrixXd::Identity(dim_dof_, dim_dof_)); // GM
        b_meas_.segment(0, dim_dof_) = Mmat_p_Ralpha_ * dq_append_.segment(3, dim_dof_ - 3);
        // Q_meas_dense.block(0, 0, dim_dof_, dim_dof_) = (Mmat_p_Ralpha_ *
        //                                                 C_GM_omega_alpha_ *
        //                                                 Mmat_p_Ralpha_.transpose())
        //                                                    .inverse();
        Q_meas_dense.block(0, 0, dim_dof_, dim_dof_) = Q_GM_meas_;
    }
    //---------------------------------------------------------------
    // LO measurement at 0
    // Q_meas_:
    //  foot_position: [   R_sb_ * (J * C_encoder_position_^{-1} * J')^{-1} * R_sb_' ];
    //  foot_velocity: [   R_sb_ * ([-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]')^{-1} * R_sb_' ];
    // b_meas_:
    //  foot_position: [  R_sb_ * fk   ];
    //  foot_velocity: [  -R_sb_ * J * dq - R_sb_ omega.cross(fk)   ];
    switch (using_lo_v_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {
        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_v_start_idx_ + i * 3,
                                                        3,
                                                        MatrixXd::Identity(3, 3)); // v

            double foot_rad = -0.023;
            // double foot_rad = 0.0;
            VectorXd p_foot = p_imu_2_foot_.block<3, 1>(i * 3, 0);
            p_foot = p_foot / p_foot.norm();
            Vector3d d = foot_rad * R_sb_ * p_foot;

            // b_meas_.segment<3>(lo_v_start_idx_ + i * 3) = -R_sb_ *
            //                                                   J_imu_2_foot_.block<3, 3>(i * 3, 0) * joint_velocity_.segment<3>(i * 3) -
            //                                               R_sb_ *
            //                                                   angular_b_.cross(p_imu_2_foot_.block<3, 1>(i * 3, 0)) +
            //                                               foot_anguler_rate_.segment<3>(i * 3).cross(d);
            b_meas_.segment<3>(lo_v_start_idx_ + i * 3) = -R_sb_ *
                                                              J_imu_2_foot_.block<3, 3>(i * 3, 0) * joint_velocity_.segment<3>(i * 3) -
                                                          R_sb_ *
                                                              angular_b_.cross(p_imu_2_foot_.block<3, 1>(i * 3, 0));
            if (contact_input_stack_.back()(i) == 0.0)
            {
                Q_meas_dense.block<3, 3>(lo_v_start_idx_ + i * 3,
                                         lo_v_start_idx_ + i * 3)
                    << Q_foot_swing_;
            }
            else
            {
                MatrixXd G_meas_i = MatrixXd::Zero(3, 9);

                MatrixXd C_meas_dense = MatrixXd::Zero(3, 3);

                G_meas_i.block<3, 3>(0, 0) = -J_imu_2_foot_.block<3, 3>(i * 3, 0);

                G_meas_i.block<3, 3>(0, 3) = -omega_skew_ * J_imu_2_foot_.block<3, 3>(i * 3, 0);

                Matrix3d kin_skew = Matrix3d::Zero();
                EigenUtils::vector3dSkew(kin_skew, p_imu_2_foot_.block<3, 1>(i * 3, 0));
                G_meas_i.block<3, 3>(0, 6) = kin_skew;

                C_meas_dense = R_sb_ * G_meas_i * C_meas_input_ *
                               G_meas_i.transpose() * R_sb_.transpose();

                Q_meas_dense.block<3, 3>(lo_v_start_idx_ + i * 3,
                                         lo_v_start_idx_ + i * 3)
                    << C_meas_dense.inverse();
            }
        }
        break;
    }
    }

    switch (using_lo_p_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {

        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        0,
                                                        -MatrixXd::Identity(3, 3)); // p
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        p_foot_start_idx_ + i * 3,
                                                        MatrixXd::Identity(3, 3)); // p_foot
            b_meas_.segment<3>(lo_p_start_idx_ + i * 3) = R_sb_ *
                                                          p_imu_2_foot_.block<3, 1>(i * 3, 0);

            Q_meas_dense.block<3, 3>(lo_p_start_idx_ + i * 3,
                                     lo_p_start_idx_ + i * 3) = R_sb_ *
                                                                (J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                                                                 C_encoder_position_ *
                                                                 J_imu_2_foot_.block<3, 3>(i * 3, 0).transpose())
                                                                    .inverse() *
                                                                R_sb_.transpose();
        }
        break;
    }
    default:
    {
        std::cout << est_type_ + " not a valid leg odom type." << std::endl;
        break;
    }
    }

    EigenUtils::SparseMatrixBlockAsignFromDense(Q_meas_, 0, 0, Q_meas_dense);

    // // A_meas_ x_i - v_{i} = b_meas_
    // mhe_qp_.addConstraints(Meas_string, b_meas_, b_meas_);
    // mhe_qp_.addConstraintDependency(Meas_string, X_T_string, A_meas_);
    // mhe_qp_.addConstraintDependency(Meas_string, V_T_string, -Identity_meas_);

    // // 0.5 * || v_{i} ||^2 _{Q_meas_}
    // mhe_qp_.addCost(Meas_string, VectorXd::Zero(dim_meas_), Q_meas_);
    // mhe_qp_.addCostDependency(Meas_string, V_T_string, Identity_meas_);

    // 0.5 * || A_meas_ x_i - b_meas_ ||^2 _{Q_meas_}
    mhe_qp_.addCost(Meas_string, b_meas_, Q_meas_);
    mhe_qp_.addCostDependency(Meas_string, X_T_string, A_meas_);

    //---------------------------------------------------------------
    // Contact Equality Constraints
    if (est_GRF_ && 1)
    {
        // dim_contact_eq_ = num_legs_ * 3;
        // A_contact_eq_.resize(dim_contact_eq_, dim_state_);
        // A_contact_eq_.setZero();

        A_contact_eq_dense_.setZero();
        b_contact_eq_ = VectorXd::Zero(dim_contact_eq_);
        // int id = 0;
        // A_contact_ineq_.setZero();
        for (int i = 0; i < num_legs_; i++)
        {
            //---------------------------------------------------------------
            // Contact Equality Constraints
            if (contact_(i) == 0.0)
            {
                // EigenUtils::SparseMatrixBlockAsignFromDense(A_contact_eq_, i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3, MatrixXd::Identity(3, 3)); // GRF
                A_contact_eq_dense_.block(i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3, 3, 3) = MatrixXd::Identity(3, 3); // GRF
                b_contact_eq_.segment<3>(i * 3) = VectorXd::Zero(3);
                // id++;
                // b_contact_ineq_lo_.segment<1>(i * 1) = -mhe_qp_.osqp_infinity * VectorXd::Ones(1);
                // b_contact_ineq_up_.segment<1>(i * 1) = mhe_qp_.osqp_infinity * VectorXd::Ones(1);
            }
            else
            {
                // EigenUtils::SparseMatrixBlockAsignFromDense(A_contact_eq_, i * 3, 9 + num_legs_ * dim_leg_ + i * 3, MatrixXd::Zero(3, 3)); // GRF

                A_contact_eq_dense_.block(i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3, 3, 3) = MatrixXd::Zero(3, 3); // GRF
                b_contact_eq_.segment<3>(i * 3) = VectorXd::Zero(3);
                // b_contact_ineq_lo_.segment<1>(i * 1) = -mhe_qp_.osqp_infinity * VectorXd::Ones(1);
                // b_contact_ineq_up_.segment<1>(i * 1) = mhe_qp_.osqp_infinity * VectorXd::Ones(1);
            }
            //---------------------------------------------------------------
            // Contact Inequality Constraints
            // EigenUtils::SparseMatrixBlockAsignFromDense(A_contact_ineq_, i * 1, 9 + num_legs_ * dim_leg_ + i * 3 + 2, MatrixXd::Identity(1, 1)); // GRF
        }
        mhe_qp_.addConstraints(Contact_string, b_contact_eq_, b_contact_eq_);
        mhe_qp_.addConstraintDependency(Contact_string, X_T_string, A_contact_eq_dense_.sparseView());

        // mhe_qp_.addConstraints(Lcp_string, b_contact_ineq_lo_, b_contact_ineq_up_);
        // mhe_qp_.addConstraintDependency(Lcp_string, X_T_string, A_contact_ineq_);
    }
    //---------------------------------------------------------------
    // Contact Complementray simplification

    // int dim_lcp = static_cast<int>(contact_input_stack_.back().sum());
    // A_lcp_.resize(dim_lcp, dim_state_);
    // A_lcp_.setZero();
    // b_lcp_.resize(dim_lcp);
    // b_lcp_.setZero();
    // int lcp_init = 0;
    // for (int i = 0; i < num_legs_; i++)
    // {
    //     if (contact_input_stack_.back()(i) == 1.0)
    //     {
    //         EigenUtils::SparseMatrixBlockAsignFromDense(A_lcp_, lcp_init, 9 + est_GRF_ * num_legs_ * (dim_leg_ + 3) + i * 3 + 2, MatrixXd::Identity(1, 1));
    //         b_lcp_.segment<1>(lcp_init) = VectorXd::Zero(1);
    //         lcp_init += 1;
    //     }
    // }
    // mhe_qp_.addConstraints(Lcp_string, b_lcp_, b_lcp_);
    // mhe_qp_.addConstraintDependency(Lcp_string, X_T_string, A_lcp_);
    // mhe_qp_.addConstraintDependency(Lcp_string, X_T_pre_string, -A_lcp_);

    mhe_qp_.updateQP(T);
}

void DecentralizedEstimation::MarginalizeKF()
{
}

// KF part validated
void DecentralizedEstimation::InitializeKF()
{
    GetMeasurement(0);

    Matrix3d R_sb = R_sb_;

    // x_prior:
    // [    0                       ;   ]   p_s
    // [    0                       ;   ]   v_s
    // [    0                       ;   ]   accel_bias_b
    // [    R_sb * fk               ;   ]   foot_position
    x_prior_.segment<3>(0) << 0.0, 0.0, 0.0;
    x_prior_.segment<3>(3) << 0.0, 0.0, 0.0;
    x_prior_.segment<3>(6) << 0.0, 0.0, 0.0;

    // C_prior_:
    //   [  C_p_init  0          0                  0;   ]
    //   [  0         C_v_init   0                  0;   ]
    //   [  0         0          C_accel_bias_init  0;   ]
    //   [  0         0          0                  C_foot_init;   ]  foot_position
    MatrixXd C_prior = MatrixXd::Zero(dim_state_, dim_state_);

    MatrixXd C_p_int = Matrix3d::Zero();
    MatrixXd C_v_int = Matrix3d::Zero();
    MatrixXd C_accel_bias_init = Matrix3d::Zero();
    MatrixXd C_foot_int = Matrix3d::Zero();
    MatrixXd C_GM_init = Matrix3d::Zero();
    MatrixXd C_GRF_init = Matrix3d::Zero();
    MatrixXd C_GM_Torso_init = MatrixXd::Zero(6, 6);

    StdVec2CovMat(params_ptr_->p_init_std_, C_p_int);
    StdVec2CovMat(params_ptr_->v_init_std_, C_v_int);
    StdVec2CovMat(params_ptr_->accel_bias_init_std_, C_accel_bias_init);
    StdVec2CovMat(params_ptr_->foot_init_std_, C_foot_int);
    StdVec2CovMat(params_ptr_->GM_init_std_, C_GM_init);
    StdVec2CovMat(params_ptr_->GRF_init_std_, C_GRF_init);
    StdVec2CovMat(params_ptr_->GM_Torso_init_std_, C_GM_Torso_init);

    C_prior.block<3, 3>(0, 0) = C_p_int;
    C_prior.block<3, 3>(3, 3) = C_v_int;
    C_prior.block<3, 3>(6, 6) = C_accel_bias_init;

    A_meas_.setZero();
    // C_meas_:
    //  foot_position: [   R_sb * J * C_encoder_position_^{-1} * J' * R_sb' ];
    //  foot_velocity: [   R_sb * [-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]' * R_sb' ];
    MatrixXd C_meas = MatrixXd::Zero(dim_meas_, dim_meas_);

    //---------------------------------------------------------------
    // GM measurement at 0
    if (est_GRF_)
    {
        Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
        Mmat_b1(Mmat_, q_append_); // To be edit to configurate euler

        VectorXd GM = VectorXd::Zero(dim_dof_);
        for (int i = 0; i < num_legs_; i++)
        {
            GM.segment(i * dim_leg_, dim_leg_) = S_.block(i * dim_leg_, 0, dim_leg_, dim_dof_) * Mmat_ * dq_append_;
        }

        // Mmat_ = [Mmat_p_  Mmat_omega_alpha]
        Mmat_p_ = Mmat_.block(0, 0, dim_dof_, 3);
        Mmat_p_Ralpha_ = Mmat_.block(0, 3, dim_dof_, dim_dof_ - 3);

        //      p   v                   accel_b_     GM_i   GRF_i
        //   [  0   - Mmat_p_        0          I       0;  ]
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    3,
                                                    -Mmat_p_); // v

        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    9,
                                                    MatrixXd::Identity(dim_dof_, dim_dof_)); // GM_rf
        b_meas_.segment(0, dim_dof_) = Mmat_p_Ralpha_ * dq_append_.segment(3, dim_dof_ - 3);
        C_meas.block(0, 0, dim_dof_, dim_dof_) = Mmat_p_Ralpha_ *
                                                 C_GM_omega_alpha_ *
                                                 Mmat_p_Ralpha_.transpose();

        x_prior_.segment(9, dim_dof_) << b_meas_.segment(0, dim_dof_);
        C_prior.block(9, 9, 6, 6) = C_GM_Torso_init;
        for (int i = 0; i < num_legs_; i++)
        {
            x_prior_.segment(9 + 6 + num_legs_ * dim_leg_ + i * 3, dim_leg_) << 0.0, 0.0, 0.0;
            C_prior.block(9 + 6 + i * dim_leg_, 9 + 6 + i * dim_leg_, dim_leg_, dim_leg_) = C_GM_init;
            C_prior.block<3, 3>(9 + 6 + num_legs_ * dim_leg_ + i * 3, 9 + 6 + num_legs_ * dim_leg_ + i * 3) = C_GRF_init;
        }
        if (est_Torso_)
        {
            MatrixXd C_Force_Torso_init = MatrixXd::Zero(6, 6);
            StdVec2CovMat(params_ptr_->Force_Torso_init_std_, C_Force_Torso_init);

            x_prior_.segment(9 + num_legs_ * (dim_leg_ + 3) + 6, 6) << VectorXd::Zero(6); // Force_Torso_prior

            C_prior.block<6, 6>(9 + 6 + num_legs_ * (dim_leg_ + 3), 9 + 6 + num_legs_ * (dim_leg_ + 3)) = C_Force_Torso_init;
        }
    }
    //---------------------------------------------------------------
    // LO measurement at 0
    // C_meas_:
    //  foot_position: [   R_sb * J * C_encoder_position_^{-1} * J' * R_sb' ];
    //  foot_velocity: [   R_sb * [-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]' * R_sb' ];

    // b_meas_:
    //  foot_position: [  R_sb * fk   ];
    //  foot_velocity: [  -R_sb * J * dq - R_sb omega.cross(fk)   ];
    switch (using_lo_v_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {
        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_v_start_idx_ + i * 3,
                                                        3,
                                                        MatrixXd::Identity(3, 3)); // v

            b_meas_.segment<3>(lo_v_start_idx_ + i * 3) = -R_sb *
                                                              J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                                                              joint_velocity_.segment<3>(i * 3) -
                                                          R_sb *
                                                              angular_b_.cross(
                                                                  p_imu_2_foot_.block<3, 1>(i * 3, 0));
            if (contact_input_stack_.back()(i) == 0.0)
            {
                C_meas.block<3, 3>(lo_v_start_idx_ + i * 3, lo_v_start_idx_ + i * 3) << C_foot_swing_;
            }
            else
            {
                MatrixXd G_meas_i = MatrixXd::Zero(3, 9);

                G_meas_i.block<3, 3>(0, 0) = -J_imu_2_foot_.block<3, 3>(i * 3, 0);

                Matrix3d omega_skew = Matrix3d::Zero();
                EigenUtils::vector3dSkew(omega_skew, angular_b_);
                G_meas_i.block<3, 3>(0, 3) = -omega_skew * J_imu_2_foot_.block<3, 3>(i * 3, 0);

                Matrix3d kin_skew = Matrix3d::Zero();
                EigenUtils::vector3dSkew(kin_skew, p_imu_2_foot_.block<3, 1>(i * 3, 0));
                G_meas_i.block<3, 3>(0, 6) = kin_skew;

                C_meas.block<3, 3>(lo_v_start_idx_ + i * 3, lo_v_start_idx_ + i * 3) = R_sb * G_meas_i * C_meas_input_ *
                                                                                       G_meas_i.transpose() * R_sb.transpose();
            }
        }
        break;
    }
    }
    switch (using_lo_p_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {

        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        0,
                                                        -MatrixXd::Identity(3, 3)); // p                                                                                             // p
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        p_foot_start_idx_ + i * 3,
                                                        MatrixXd::Identity(3, 3)); // p_foot

            b_meas_.segment<3>(lo_p_start_idx_ + i * 3) = R_sb *
                                                          p_imu_2_foot_.block<3, 1>(i * 3, 0);
            C_meas.block<3, 3>(lo_p_start_idx_ + i * 3,
                               lo_p_start_idx_ + i * 3) = R_sb *
                                                          J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                                                          C_encoder_position_ *
                                                          J_imu_2_foot_.block<3, 3>(i * 3, 0).transpose() *
                                                          R_sb.transpose();
            C_prior.block<3, 3>(p_foot_start_idx_ + i * 3, p_foot_start_idx_ + i * 3) = C_foot_int;
        }
        x_prior_.segment(p_foot_start_idx_, 3 * num_legs_) = b_meas_.segment(lo_p_start_idx_,
                                                                             3 * num_legs_);
        break;
    }
    default:
    {
        break;
    }
    }

    // KF prior setup at time 0
    x_KF_ = x_prior_;
    C_KF_ = C_prior;

    // KF measurement correction at time 0
    K_KF_ = C_KF_ * A_meas_.transpose() * (A_meas_ * C_KF_ * A_meas_.transpose() + C_meas).inverse();
    x_KF_ = x_KF_ + K_KF_ * (b_meas_ - A_meas_ * x_KF_);
    C_KF_ = (MatrixXd::Identity(dim_state_, dim_state_) - K_KF_ * A_meas_) * C_KF_;
}

void DecentralizedEstimation::UpdateKF()
{
    //---------------------------------------------------------------
    // KF prediction update

    Matrix3d R_sb = R_sb_;
    Vector3d accel_s = accel_s_;
    double dt = dt_;

    // b_dyn_:
    // [    - 0.5 * dt^2 * accel_s; ]   p_s
    // [    - dt * accel_s        ; ]   v_s
    // [    0                     ; ]   accel_bias_b
    // [    0                     ; ]   p_foot_s
    b_dyn_.segment(0, 3) << -0.5 * dt * dt * accel_s;
    b_dyn_.segment(3, 3) << -dt * accel_s;

    // A_dyn_:
    // p_{i+1} = p_i + dt * v_i + 0.5 * dt^2 * accel_s - 0.5 * dt^2 * R_sb * accel_b_bias + w_p;
    // v_{i+1} = v_i + dt * accel_s - dt * R_sb * accel_b_bias + w_v;
    // accel_b_bias_{i+1} = accel_b_bias_{i} + w_a_bias
    // p_foot_{i+1} = p_foot_{i} + w_foot
    // [    I   dt* I   - 0.5 * dt^2 * R_sb;  ]
    // [    0   I       - dt * R_sb        ;  ]
    // [    0   0       I                  ;  ]
    A_dyn_.setIdentity();
    EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_, 0, 3, dt * Matrix3d::Identity());
    EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_, 0, 6, -dt * dt / 2 * R_sb);
    EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_, 3, 6, -dt * R_sb);

    MatrixXd C_dyn = MatrixXd::Zero(dim_state_, dim_state_);

    // G_dyn:
    // [  R_sb * dt          -0.5 * R_sb *dt^2   0           0;      ]
    // [  0                  -R_sb *dt           0           0;      ]
    // [  0                  0                   I * dt      0;      ]
    // [  0                  0                   0           R_sb * dt;]
    // C_dyn = G_dyn * diag[C_velocity, C_accel_, C_accel_bias_, C_foot_  ] * G_dyn'
    // C_foot = infinite_, if contact(feet) = false
    MatrixXd G_dyn = MatrixXd::Zero(dim_state_, dim_state_);
    G_dyn.block<3, 3>(0, 0) = R_sb * dt;
    G_dyn.block<3, 3>(0, 3) = -0.5 * R_sb * dt * dt;
    G_dyn.block<3, 3>(3, 3) = -R_sb * dt;
    G_dyn.block<3, 3>(6, 6) = Matrix3d::Identity() * dt;

    MatrixXd C_input = MatrixXd::Zero(dim_state_, dim_state_);
    C_input.block<3, 3>(0, 0) = C_p_;
    C_input.block<3, 3>(3, 3) = C_accel_;
    C_input.block<3, 3>(6, 6) = C_accel_bias_;

    //---------------------------------------------------------------
    // GM, GRF states update
    // GM_i^+ = GM_i + dt* S_i_ * C^T_p * v + dt * S_i_ * J_i * GRF_i + dt * {S_i_ * C^T_{w,alpha} - S_i_ * G  + B * u}
    // GRF_i^+ = GRF_i
    if (est_GRF_)
    {
        Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
        Mmat_b1(Mmat_, q_append_); // To be edit to configurate euler

        Cmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);

        dq_append_.segment<3>(0) = x_KF_.segment<3>(3);

        SymFunction::C_mat_go1(Cmat_, q_append_, dq_append_);
        Cmat_T_p_ = (Cmat_.block(0, 0, 3, dim_dof_)).transpose();

        Bmat_ = MatrixXd::Zero(dim_dof_, dim_leg_ * num_legs_);
        SymFunction::u_map_b1_description(Bmat_, q_append_);

        Gvec_ = MatrixXd::Zero(dim_dof_, 1);
        SymFunction::Ge_vec_b1_description(Gvec_, q_append_, VectorXd::Zero(dim_dof_));
        Gvec_ = -Gvec_;

        int GM_start_idx = 9;
        int GRF_foot_start_idx = 9 + 6 + num_legs_ * dim_leg_;

        // A_dyn_:
        //      p   v                 accel_b_     GM_i     GRF_i
        //   [  0   dt *  C^T_p       0            I        dt  * J_i^T;  ]
        EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_,
                                                    GM_start_idx,
                                                    3,
                                                    dt * Cmat_T_p_); // v

        b_dyn_.segment(GM_start_idx, dim_dof_) = -dt *
                                                 (Cmat_.block(3, 0, dim_dof_ - 3, dim_dof_).transpose() *
                                                      dq_append_.segment(3, dim_dof_ - 3) -
                                                  Gvec_ +
                                                  Bmat_ * joint_effort_);
        G_dyn.block(GM_start_idx, GM_start_idx, 6, 6) = MatrixXd::Identity(6, 6) * dt;
        C_input.block(GM_start_idx,
                      GM_start_idx,
                      6,
                      6) = C_GM_Torso_process_;

        for (int i = 0; i < num_legs_; i++)
        {

            MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
            JFunctions[i](J_i, q_append_);
            EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_,
                                                        GM_start_idx,
                                                        GRF_foot_start_idx + i * 3,
                                                        dt * J_i.transpose());
            // b_dyn_:
            // - dt * {S_i_ * C^T_{pitch,alpha} [w;dalpha] - S_i_ * G  + S_i_ * B * u}

            G_dyn.block(GM_start_idx + 6 + i * dim_leg_,
                        GM_start_idx + 6 + i * dim_leg_,
                        dim_leg_,
                        dim_leg_) = MatrixXd::Identity(dim_leg_, dim_leg_) * dt;
            G_dyn.block(GRF_foot_start_idx + i * 3, GRF_foot_start_idx + i * 3, 3, 3) = MatrixXd::Identity(3, 3) * dt;
            C_input.block(GM_start_idx + 6 + i * dim_leg_,
                          GM_start_idx + 6 + i * dim_leg_,
                          dim_leg_,
                          dim_leg_) = C_GM_process_;
            C_input.block(GRF_foot_start_idx + i * 3,
                          GRF_foot_start_idx + i * 3,
                          3,
                          3) = C_GRF_process_;
        }

        if (est_Torso_)
        {
            MatrixXd J_torso = MatrixXd::Zero(6, dim_dof_);
            J_torso.block(0, 0, 6, 6) = MatrixXd::Identity(6, 6);
            EigenUtils::SparseMatrixBlockAsignFromDense(A_dyn_,
                                                        GM_start_idx, wrench_torso_start_idx_,
                                                        dt * J_torso.transpose()); // Wrench
            G_dyn.block(wrench_torso_start_idx_, wrench_torso_start_idx_, 6, 6) = MatrixXd::Identity(6, 6) * dt;
            C_input.block(wrench_torso_start_idx_, wrench_torso_start_idx_, 6, 6) = C_Force_Torso_process_;
        }
    }
    //---------------------------------------------------------------
    // Leg states update
    switch (using_lo_p_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {
        for (int i = 0; i < num_legs_; i++)
        {
            G_dyn.block<3, 3>(p_foot_start_idx_ + i * 3, p_foot_start_idx_ + i * 3) = R_sb * dt;

            // left foot covariance
            if (contact_input_stack_.back()(i) == 0.0) // if contact at foot_idx, fixed the stance foot position using small C_foot_slide
            {
                C_input.block<3, 3>(p_foot_start_idx_ + 3 * i, p_foot_start_idx_ + 3 * i) << C_foot_swing_;
            }
            else
            {
                C_input.block<3, 3>(p_foot_start_idx_ + 3 * i, p_foot_start_idx_ + 3 * i) = C_foot_slide_;
            }
        }
        break;
    }
    default:
    {
        break;
    }
    }

    GetMeasurement(0);

    //---------------------------------------------------------------
    // Impact correction
    if (est_GRF_ && using_impact_map_)
    {
        int impact_count = (contact_change_stack_.back().array() == 1.0).count();
        MatrixXd J_impact_B = MatrixXd::Zero(dim_dof_, 3 * impact_count);
        MatrixXd J_impact_C = MatrixXd::Zero(3 * impact_count, dim_dof_);
        int foot_cout = 0;
        MatrixXd Mmat_inv = Mmat_.inverse();

        for (int i = 0; i < num_legs_; i++)
        {
            if (contact_change_stack_.back()(i) == 1.0)
            {
                MatrixXd J_i = MatrixXd::Zero(3, dim_dof_);
                JFunctions[i](J_i, q_append_);
                J_impact_B.block(0, 0 + foot_cout * 3, dim_dof_, 3) = -J_i.transpose();
                J_impact_C.block(0 + foot_cout * 3, 0, 3, dim_dof_) = J_i * Mmat_inv;
                foot_cout += 1;
            }
        }

        MatrixXd A_impact = MatrixXd::Identity(dim_dof_, dim_dof_) - J_impact_B * (J_impact_C * J_impact_B).inverse() * J_impact_C;
        std::cout << A_impact << std::endl;
        SparseMatrix<double> A_impact_reorder;
        A_impact_reorder.resize(dim_state_, dim_state_);
        A_impact_reorder.setIdentity();
        EigenUtils::SparseMatrixBlockReplaceFromDense(A_impact_reorder,
                                                      9,
                                                      9,
                                                      A_impact);
        A_dyn_ = A_dyn_ * A_impact_reorder;

        // MatrixXd F_map = (J_impact_C * J_impact_B).inverse() * J_impact_C;
        // VectorXd dq_local = dq_append_;
        // dq_local.segment<3>(0) = x_KF_.segment<3>(3);
        // std::cout << F_map * Mmat_ * dq_local << std::endl;
    }

    x_KF_ = A_dyn_ * x_KF_ - b_dyn_;
    C_dyn = G_dyn * C_input * G_dyn.transpose();
    C_KF_ = A_dyn_ * C_KF_ * A_dyn_.transpose() + C_dyn;
    //---------------------------------------------------------------
    // KF measurement correction
    R_sb = R_sb_;

    // C_meas_:
    //  foot_position: [   R_sb * J * C_encoder_position_^{-1} * J' * R_sb' ];
    //  foot_velocity: [   R_sb * [-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]' * R_sb' ];
    MatrixXd C_meas = MatrixXd::Zero(dim_meas_, dim_meas_);
    A_meas_.setZero();

    //---------------------------------------------------------------
    // GM measurement at T
    if (est_GRF_)
    {
        Mmat_ = MatrixXd::Zero(dim_dof_, dim_dof_);
        Mmat_b1(Mmat_, q_append_); // To be edit to configurate euler

        VectorXd GM = VectorXd::Zero(dim_dof_);
        for (int i = 0; i < num_legs_; i++)
        {
            GM.segment(i * dim_leg_, dim_leg_) = S_.block(i * dim_leg_, 0, dim_leg_, dim_dof_) * Mmat_ * dq_append_;
        }

        // Mmat_ = [Mmat_p_  Mmat_omega_alpha]
        Mmat_p_ = Mmat_.block(0, 0, dim_dof_, 3);
        Mmat_p_Ralpha_ = Mmat_.block(0, 3, dim_dof_, dim_dof_ - 3);

        //      p   v                   accel_b_     GM_i   GRF_i
        //   [  0   - Mmat_p_        0          I       0;  ]
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    3,
                                                    -Mmat_p_); // v
        EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                    0,
                                                    9,
                                                    MatrixXd::Identity(dim_dof_, dim_dof_)); // GM_rf
        b_meas_.segment(0, dim_dof_) = Mmat_p_Ralpha_ * dq_append_.segment(3, dim_dof_ - 3);
        C_meas.block(0, 0, dim_dof_, dim_dof_) = Mmat_p_Ralpha_ *
                                                 C_GM_omega_alpha_ *
                                                 Mmat_p_Ralpha_.transpose();
    }

    //---------------------------------------------------------------
    // LO measurement at 0
    // C_meas_:
    //  foot_position: [   R_sb * J * C_encoder_position_^{-1} * J' * R_sb' ];
    //  foot_velocity: [   R_sb * [-J, -omega^x * J, fk^x] *
    //                          [C_encoder_position_,C_encoder_velocity_,C_gyro_] *
    //                          [-J, -omega^x * J, fk^x]' * R_sb' ];

    // b_meas_:
    //  foot_position: [  R_sb * fk   ];
    //  foot_velocity: [  -R_sb * J * dq - R_sb omega.cross(fk)   ];
    switch (using_lo_v_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {
        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_, lo_v_start_idx_ + i * 3, 3, MatrixXd::Identity(3, 3)); // v
            double foot_rad = -0.023;
            // double foot_rad = 0.0;
            VectorXd p_foot = p_imu_2_foot_.block<3, 1>(i * 3, 0);
            p_foot = p_foot / p_foot.norm();
            Vector3d d = foot_rad * R_sb * p_foot;

            // std::cout << foot_anguler_rate_.segment<3>(i * 3).cross(d) << std::endl;
            // b_meas_.segment<3>(lo_v_start_idx_ + i * 3) = -R_sb *
            //                                                   J_imu_2_foot_.block<3, 3>(i * 3, 0) *
            //                                                   joint_velocity_stack_.back().segment<3>(i * 3) -
            //                                               R_sb *
            //                                                   angular_b_input_stack_.back().cross(
            //                                                       p_imu_2_foot_.block<3, 1>(i * 3, 0)) +
            //                                               foot_anguler_rate_.segment<3>(i * 3).cross(d);
            b_meas_.segment<3>(lo_v_start_idx_ + i * 3) = -R_sb *
                                                              J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                                                              joint_velocity_.segment<3>(i * 3) -
                                                          R_sb *
                                                              angular_b_.cross(
                                                                  p_imu_2_foot_.block<3, 1>(i * 3, 0));
            if (contact_input_stack_.back()(i) == 0.0)
            {
                C_meas.block<3, 3>(lo_v_start_idx_ + i * 3, lo_v_start_idx_ + i * 3) << C_foot_swing_;
            }
            else
            {
                MatrixXd G_meas_i = MatrixXd::Zero(3, 9);

                G_meas_i.block<3, 3>(0, 0) = -J_imu_2_foot_.block<3, 3>(i * 3, 0);

                Matrix3d omega_skew = Matrix3d::Zero();
                EigenUtils::vector3dSkew(omega_skew, angular_b_);
                G_meas_i.block<3, 3>(0, 3) = -omega_skew * J_imu_2_foot_.block<3, 3>(i * 3, 0);

                Matrix3d kin_skew = Matrix3d::Zero();
                EigenUtils::vector3dSkew(kin_skew, p_imu_2_foot_.block<3, 1>(i * 3, 0));
                G_meas_i.block<3, 3>(0, 6) = kin_skew;

                C_meas.block<3, 3>(lo_v_start_idx_ + i * 3, lo_v_start_idx_ + i * 3) = R_sb * G_meas_i * C_meas_input_ *
                                                                                       G_meas_i.transpose() * R_sb.transpose();
            }
        }
        break;
    }
    }

    switch (using_lo_p_)
    {
    case 0:
    {
        break;
    }
    case 1:
    {
        for (int i = 0; i < num_legs_; i++)
        {
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        0,
                                                        -MatrixXd::Identity(3, 3)); // p                                                                                             // p
            EigenUtils::SparseMatrixBlockAsignFromDense(A_meas_,
                                                        lo_p_start_idx_ + i * 3,
                                                        p_foot_start_idx_ + i * 3,
                                                        MatrixXd::Identity(3, 3)); // p_foot

            b_meas_.segment<3>(lo_p_start_idx_ + i * 3) = R_sb *
                                                          p_imu_2_foot_.block<3, 1>(i * 3, 0);
            C_meas.block<3, 3>(lo_p_start_idx_ + i * 3,
                               lo_p_start_idx_ + i * 3) = R_sb *
                                                          J_imu_2_foot_.block<3, 3>(i * 3, 0) *
                                                          C_encoder_position_ *
                                                          J_imu_2_foot_.block<3, 3>(i * 3, 0).transpose() *
                                                          R_sb.transpose();
        }
        break;
    }
    default:
    {
        break;
    }
    }

    // KF measurement correction
    K_KF_ = C_KF_ * A_meas_.transpose() * (A_meas_ * C_KF_ * A_meas_.transpose() + C_meas).inverse();
    x_KF_ = x_KF_ + K_KF_ * (b_meas_ - A_meas_ * x_KF_);
    C_KF_ = (MatrixXd::Identity(dim_state_, dim_state_) - K_KF_ * A_meas_) * C_KF_;
}

// Tested
void DecentralizedEstimation::GetMeasurement(int T)
{

    foot_anguler_rate_ = robot_sub_ptr_->foot_anguler_rate_;

    Quaterniond q_bi;
    // q_bi.w() = 0.9998;
    // q_bi.x() = -0.020;
    // q_bi.y() = -0.0035;
    // q_bi.z() = -0.0192;

    q_bi.w() = 1;
    q_bi.x() = 0;
    q_bi.y() = 0;
    q_bi.z() = 0;

    Matrix3d R_bi = q_bi.normalized().toRotationMatrix();

    // get current measurement from robot_sub_ptr_
    R_sb_ = (robot_sub_ptr_->quaternion_.normalized()).toRotationMatrix();
    EigenUtils::QuaternionToEuler(robot_sub_ptr_->quaternion_.normalized(), euler_sb_);

    double imu_time_ = robot_sub_ptr_->imu_time_;
    accel_b_ = R_bi * robot_sub_ptr_->accel_b_; // returns the accelerometer readings in body frame;
    accel_s_ = R_sb_ * accel_b_ + gravity_;
    angular_b_ = R_bi * robot_sub_ptr_->angular_b_;

    EigenUtils::vector3dSkew(omega_skew_, angular_b_);

    EigenUtils::AngularToEulerRates(angular_b_, euler_sb_, euler_rate_sb_);

    joint_position_ = robot_sub_ptr_->joint_states_position_;
    joint_velocity_ = robot_sub_ptr_->joint_states_velocity_;
    joint_effort_ = robot_sub_ptr_->joint_states_effort_;

    p_imu_2_foot_ = robot_sub_ptr_->p_imu_2_foot_;
    J_imu_2_foot_ = robot_sub_ptr_->J_imu_2_foot_;
    contact_ = robot_sub_ptr_->contact_;

    q_append_.segment<3>(3) = euler_sb_;
    q_append_.segment(6, dim_leg_ * num_legs_) = joint_position_;

    dq_append_.setZero();
    dq_append_.segment<3>(3) = euler_rate_sb_;
    dq_append_.segment(6, dim_leg_ * num_legs_) = joint_velocity_;

    if (contact_input_stack_.size() != 0)
    {
        contact_change_ = contact_ - contact_input_stack_.back();
    }
    else
    {
        contact_change_ = VectorXd::Zero(num_legs_);
    }
    //------------------------------------------------------
    // Sychronize the VO frames to the nearest IMU frames to the left
    vo_new_meas_flag_ = robot_sub_ptr_->vo_new_;
    if (vo_new_meas_flag_ && imu_time_stack_.size() > 0) // log when new vo transformation was subscribed
    {

        vo_p_body_pre_2_body_ = robot_sub_ptr_->vo_p_body_pre_2_body_; // Relative translation between body frame in body_pre

        double vo_time_pre_ = robot_sub_ptr_->vo_time_pre_;
        double vo_time_now_ = robot_sub_ptr_->vo_time_now_;
        robot_sub_ptr_->vo_new_ = false; // meassage logged

        // sychronize VO timestamp to the first IMU timestamp to the left
        // IMU time is the rclcpp::clock.now() when the imu_callback gets callled in subscriber callback
        auto idx_pre_ptr = std::upper_bound(imu_time_stack_.begin(),
                                            imu_time_stack_.end(), vo_time_pre_); // find the first imu time bigger than vo_time_pre_

        if (idx_pre_ptr == imu_time_stack_.begin())
        {
            std::cout << "not storing enough imu info, failed to interpolate" << std::endl;
            // could happen at the beginning of the MHE
            // discard the vo_meas
            // this if condition should only happens at during the initalization, so give some time for the imu to store enough info
        }
        else
        {
            int imu_sychron_idx_pre = std::distance(imu_time_stack_.begin(), idx_pre_ptr) - 1;
            // double imu_sychron_time_pre = imu_time_stack_[imu_sychron_idx_pre];
            R_vo_sb_pre_ = R_input_rotation_stack_[imu_sychron_idx_pre]; // R_world_2_vo_pre, used to get relative translation between vo frames in world coordinates

            auto idx_now_ptr = std::upper_bound(imu_time_stack_.begin(),
                                                imu_time_stack_.end(), vo_time_now_); // find the first imu time bigger than vo_time_
            int imu_sychron_idx_now = std::distance(imu_time_stack_.begin(), idx_now_ptr) - 1;

            p_vo_accmulate_ += R_vo_sb_pre_ * vo_p_body_pre_2_body_; // accumulate translation in world frame, (exist a constant offset)

            int imu_window_start_idx = int(imu_time_stack_.size() - std::min(N_, T));
            int imu_interpolation_start_idx = std::max(imu_window_start_idx, imu_sychron_idx_pre); // interpolation starts from either the start of the window or image_frame_sychron_pre
            double time_interpolate_start = imu_time_stack_[imu_interpolation_start_idx];
            int discrete_time_interpolate_start = discrete_time_stack[imu_interpolation_start_idx]; // the discrete time used to find the appropriate VO_measurement_discrete_time

            // assuming the control points are uniform in time, the assumption is reasonable when using RealSense D455 + ORBSLAM3
            vo_curve_.add_way_point(p_vo_accmulate_, vo_time_now_); // the control point is consists of P_now, P_pre, P_{pre-1}, P_{pre-2}

            if (imu_sychron_idx_now > imu_window_start_idx && vo_curve_._way_points.size() >= 4)
            {
                double insert_relative_idx = imu_interpolation_start_idx - imu_window_start_idx; // the realtive idx in the window/horizon that the imterpolation starts
                double interpolate_num = imu_sychron_idx_now - imu_interpolation_start_idx + 1;  // interpolate using a constant time difference 1/N, also approximation

                vo_curve_.set_interval(time_interpolate_start, interpolate_num, dt_);

                // need to decide the knot vector here, otherwise the control point is treated as equaly spaced in time
                vo_curve_.interpolate_waypoint();

                vo_insert_idx_stack_.push_back(insert_relative_idx);
                vo_insert_discrete_time_stack_.push_back(discrete_time_interpolate_start);

                vo_to_be_processed_flag_ = true;
            }
            // vo_pose_body_pre_2_body_stack.push_back(p_vo_accmulate_);
            // R_vo_sb_rotation_stack.push_back(R_vo_sb_pre_); // record the R_world_2_body
            // vo_sychron_time_stack_.push_back(imu_sychron_time_pre);
            // vo_time_stack_.push_back(vo_time_pre_);
        }
    }

    //------------------------------------------------------
    // Stack and un-Stack
    imu_time_stack_.push_back(imu_time_);
    discrete_time_stack.push_back(T);
    R_input_rotation_stack_.push_back(R_sb_);
    // euler_input_rotation_stack_.push_back(euler_sb_);
    // angular_b_input_stack_.push_back(angular_b_);
    // accel_s_input_stack_.push_back(accel_s_);

    // joint_position_stack_.push_back(joint_position_);
    // joint_velocity_stack_.push_back(joint_velocity_);
    // joint_effort_stack_.push_back(joint_effort_);

    // p_imu_2_foot_stack_.push_back(p_imu_2_foot_);
    // J_imu_2_foot_stack_.push_back(J_imu_2_foot_);
    contact_input_stack_.push_back(contact_);
    contact_change_stack_.push_back(contact_change_);

    gt_v_s_ = robot_sub_ptr_->gt_v_s_;
    // gt_v_s_input_stack_.push_back(gt_v_s_);

    // Only keep the previous N_ measurements, but the storage size is arbitary since no stack.front is used
    // discrete time now: T;
    // discrete time optimization window start: T-N_
    // if (int(accel_s_input_stack_.size()) > N_ + 1) // should be N_+1 to store the measments in the current optimization window
    if (int(accel_s_input_stack_.size()) > 4 * N_ + 1)
    {
        imu_time_stack_.erase(imu_time_stack_.begin());
        discrete_time_stack.erase(discrete_time_stack.begin());
        R_input_rotation_stack_.erase(R_input_rotation_stack_.begin());
        // euler_input_rotation_stack_.erase(euler_input_rotation_stack_.begin());
        // angular_b_input_stack_.erase(angular_b_input_stack_.begin());
        // accel_s_input_stack_.erase(accel_s_input_stack_.begin());

        // joint_position_stack_.erase(joint_position_stack_.begin());
        // joint_velocity_stack_.erase(joint_velocity_stack_.begin());
        // joint_effort_stack_.erase(joint_effort_stack_.begin());

        // p_imu_2_foot_stack_.erase(p_imu_2_foot_stack_.begin());
        // J_imu_2_foot_stack_.erase(J_imu_2_foot_stack_.begin());
        contact_input_stack_.erase(contact_input_stack_.begin());
        contact_change_stack_.erase(contact_change_stack_.begin());

        // gt_v_s_input_stack_.erase(gt_v_s_input_stack_.begin());
    }

    if (int(vo_insert_idx_stack_.size()) >= N_ + 1)
    {
        // R_vo_sb_rotation_stack_.erase(R_vo_sb_rotation_stack.begin());
        // vo_sychron_time_stack_.erase(vo_sychron_time_stack_.begin());
        // vo_time_stack_.erase(vo_time_stack_.begin());
        vo_insert_idx_stack_.erase(vo_insert_idx_stack_.begin());
        vo_insert_discrete_time_stack_.erase(vo_insert_discrete_time_stack_.begin());
    }
}

void DecentralizedEstimation::UpdateVOConstraints(int T)
{
    std::map<int, Vector3d> vo_constraints_idx_regs;
    std::map<int, double> vo_reliability_idx_regs;
    double vo_reliability = 1.0;

    for (int i = 0; i < vo_curve_.node_count() - 1; ++i)
    {
        Vector3d pose_world_idx = vo_curve_._distances[i + 1];

        // sparse meas
        // int idx = (vo_insert_idx_stack_.back() + i) *
        //               (dim_meas_ + dim_state_ + dim_cam_) +
        //           dim_meas_ + dim_state_;

        // dense meas
        int idx = (vo_insert_idx_stack_.back() + i) *
                      (dim_state_ + dim_cam_) +
                  dim_state_;
        vo_constraints_idx_regs.insert({idx, pose_world_idx});
        vo_reliability_idx_regs.insert({vo_insert_discrete_time_stack_.back() + i, vo_reliability});

        // update the constraintbound for the marginalization, no marginalization if unbounded
        std::string update_name = "VO_measurement_" + std::to_string(vo_insert_discrete_time_stack_.back() + i);
        mhe_qp_.updateConstraintBound(update_name, -pose_world_idx, -pose_world_idx, true);
        // mhe_qp_.updateCostGain(update_name, vo_reliability);
    }
    mhe_qp_.Update_Image_bound(vo_constraints_idx_regs, vo_reliability_idx_regs);
}

void DecentralizedEstimation::reset()
{

    mhe_qp_.resetQP();
}

void DecentralizedEstimation::Mmat_b1(MatrixXd &Mmat, const VectorXd &q)
{
    assert_size_matrix(q, dim_dof_, 1);
    assert_size_matrix(Mmat, dim_dof_, dim_dof_);

    Mmat = MatrixXd::Zero(dim_dof_, dim_dof_);

    // for (int i = 0; i < numFunctions; ++i)
    for (int i = 0; i < 14; ++i)
    {
        MatrixXd Mmat_temp = MatrixXd::Zero(dim_dof_, dim_dof_);
        MmatFunctions[i](Mmat_temp, q);
        Mmat += Mmat_temp;
    }
}

void DecentralizedEstimation::StdVec2GainMat(const std::vector<double> &std, MatrixXd &Gain)
{
    // Ensure the matrix is square and its size matches the size of the std vector
    assert(Gain.rows() == Gain.cols());
    assert(Gain.rows() == std.size());

    Gain = MatrixXd::Zero(Gain.rows(), Gain.cols());

    for (size_t i = 0; i < std.size(); ++i)
    {
        Gain(i, i) = 1 / std::pow(std[i], 2);
    }
}

void DecentralizedEstimation::StdVec2CovMat(const std::vector<double> &std, MatrixXd &Cov)
{
    // Ensure the matrix is square and its size matches the size of the std vector
    assert(Cov.rows() == Cov.cols());
    assert(Cov.rows() == std.size());

    Cov = MatrixXd::Zero(Cov.rows(), Cov.cols());

    for (size_t i = 0; i < std.size(); ++i)
    {
        Cov(i, i) = std::pow(std[i], 2);
    }
}

void DecentralizedEstimation::tic(std::string str, int mode)
{
    static std::chrono::_V2::system_clock::time_point t_start;

    if (mode == 0)
        t_start = std::chrono::high_resolution_clock::now();
    else
    {

        auto t_end = std::chrono::high_resolution_clock::now();
        solve_time += (t_end - t_start).count() * 1E-9;

        std::cout << str + " elapsed time: " << solve_time / solve_num << " seconds\n";
        solve_num += 1;
    }
}

void DecentralizedEstimation::toc(std::string str) { tic(str, 1); }
