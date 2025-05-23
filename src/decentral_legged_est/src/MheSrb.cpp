#include "decentral_legged_est/MheSrb.hpp"
#include <fstream>
#include <thread>
#include <chrono>
using namespace MHEutils;

#include <iostream>
#include <fstream>

// define the format you want, you only need one instance of this...
// see https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

// writing functions taking Eigen types as parameters,
// see https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
template <typename Derived>
void writeToCSVfile(std::string name, const Eigen::MatrixBase<Derived> &matrix)
{
    std::ofstream file(name.c_str(), std::ios::trunc); // trunc mode ensures file is cleared

    if (file.is_open())
    {
        // Define CSV format
        file << matrix.format(CSVFormat);
        std::cout << "File written successfully: " << name << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << name << std::endl;
    }
}

const static Eigen::IOFormat CSVFormat1(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

// Writing function for Eigen::VectorXd
template <typename Derived>
void writeVectorToCSVfile(std::string name, const Eigen::MatrixBase<Derived> &vector)
{
    // Open the file in truncate mode to clear previous contents
    std::ofstream file(name.c_str(), std::ios::trunc); // trunc mode ensures the file is cleared

    if (file.is_open())
    {
        // Write the vector in CSV format
        file << vector.format(CSVFormat1);
        std::cout << "File written successfully: " << name << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file: " << name << std::endl;
    }
}

MHEproblem::MHEproblem()
{
}

void MHEproblem::setHorizon(int N, int dim_state, int dim_meas, int dim_cam)
{
    this->N = N;

    this->dim_state = dim_state;
    this->dim_meas = dim_meas;
    this->dim_cam = dim_cam;
    // this->dim_contact = dim_contact;

    I.resize(dim_state, dim_state);
    I.setIdentity();

    I_meas.resize(dim_meas, dim_meas);
    I_meas.setIdentity();

    I_dyn.resize(dim_state, dim_state);
    I_dyn.setIdentity();

    I_corelate.resize(dim_state + dim_cam, dim_state + dim_cam);
    I_corelate.setIdentity();

    I_mag.resize(dim_state + dim_cam + dim_meas, dim_state + dim_cam + dim_meas);
    I_mag.setIdentity();

    g = VectorXd::Zero(dim_state);
}

void MHEproblem::addVariable(std::string VarName_, int VarSize_)
{
    // var0 is the initial guess
    var_idx_regs.insert({VarName_, nVarEnd}); // add the variable and the idx where it starts (nVarEnd) in the opt variable X.
    nVarEnd += VarSize_;                      // nVarEnd is the ending index of current interval
    nVar = nVarEnd - nVarStart;               // nVar is the current length of the interval
}

// 0.5 * ||Ax - b||_Q
// Add the cost term to Cost_new_regs
void MHEproblem::addCost(std::string CostName_, VectorXd b, SparseMatrix<double> Q)
{
    nCost += b.size();
    nCost_new += b.size();
    struct CostItem temp;
    temp.CostName = CostName_;
    temp.b = b;
    temp.Q = Q;
    Cost_new_regs.insert({CostName_, temp});
}

// lb <= Ax <= ub
// Add constraints to LinearConstraint_new_regs & record the size
void MHEproblem::addConstraints(std::string ConstraintName_, VectorXd lb, VectorXd ub)
{
    nConstraints += ub.size();
    nConstraints_new += ub.size();
    struct LinearConstraint temp;
    temp.CostraintName = ConstraintName_;
    temp.lb = lb;
    temp.ub = ub;
    LinearConstraint_new_regs.insert({ConstraintName_, temp});
    LinearConstraint_new_name_vec.push_back(ConstraintName_); // used to assemble the constraints by the order they are added
}

void MHEproblem::addTimeStep(int Time, int Var_Size, int Constraint_Size)
{

    struct TimeStep temp;

    temp.Time = Time;
    temp.Var_Size = Var_Size;
    temp.Constraint_Size = Constraint_Size;
    Time_regs.insert({Time, temp});
}

// Assemble constraints in LinearConstraint_regs
void MHEproblem::assembleConstraints()
{
    A_constraint_Matrices.clear();
    lb_vectors.clear();
    ub_vectors.clear();

    for (const auto &eachcon : LinearConstraint_regs)
    {

        MatrixXd Aconst = MatrixXd::Zero(eachcon.second.lb.size(), nVar);

        for (const auto &const_Ax : eachcon.second.depVarMap)
        {
            std::string var_name = const_Ax.first;
            // sparse to dense
            MatrixXd A = MatrixXd(const_Ax.second);
            if (var_idx_regs.count(var_name))
            {
                Aconst.block(0, var_idx_regs[var_name] - nVarStart, A.rows(), A.cols()) << A;
            }
        }
        A_constraint_Matrices.push_back(Aconst.sparseView());
        lb_vectors.push_back(eachcon.second.lb);
        ub_vectors.push_back(eachcon.second.ub);
        Aconst.setZero();
    }

    // Note that the matrix is sparse
    Aconstrsparse.resize(nConstraints, nVar);
    Aconstrsparse.setZero();

    // constraint bounds;
    lb_all.resize(nConstraints);
    ub_all.resize(nConstraints);
    lb_all.setZero();
    ub_all.setZero();

    int rowIdx = 0;
    for (const auto &Ai : A_constraint_Matrices)
    {
        // Aconstrsparse.block(rowIdx, 0, Ai.rows(), Ai.cols()) << Ai;
        EigenUtils::SparseMatrixBlockAsign(Aconstrsparse, rowIdx, 0, Ai);
        rowIdx += Ai.rows();
    }

    rowIdx = 0;
    for (const auto &lb_i : lb_vectors)
    {
        lb_all.segment(rowIdx, lb_i.size()) << lb_i;
        rowIdx += lb_i.rows();
    }

    rowIdx = 0;
    for (const auto &ub_i : ub_vectors)
    {
        ub_all.segment(rowIdx, ub_i.size()) << ub_i;
        rowIdx += ub_i.rows();
    }
}

// Assemble costs in Cost_regs
void MHEproblem::assembleCost()
{
    A_cost_Matrices.clear();
    b_cost_Vectors.clear();
    Q_cost_Matrices.clear();

    for (const auto &eachcost : Cost_regs)
    {

        MatrixXd Acost = MatrixXd::Zero(eachcost.second.b.size(), nVar);

        for (const auto &cost_Ax : eachcost.second.depVarMap)
        {
            std::string var_name = cost_Ax.first;
            MatrixXd A = MatrixXd(cost_Ax.second);
            if (var_idx_regs.count(var_name))
            {
                Acost.block(0, var_idx_regs[var_name] - nVarStart, A.rows(), A.cols()) << A;
            }
        }
        A_cost_Matrices.push_back(Acost.sparseView());
        b_cost_Vectors.push_back(eachcost.second.b);
        Q_cost_Matrices.push_back(eachcost.second.Q);
        Acost.setZero();
    }

    // Note that the matrix is sparse
    A_cost_all.resize(nCost, nVar); // problem with eigen on sparse matrix, which cannot do block assignment
    A_cost_all.setZero();

    b_cost_all.resize(nCost);
    b_cost_all.setZero();

    Q_cost_all.resize(nCost, nCost); // problem with eigen on sparse matrix, which cannot do block assignment
    Q_cost_all.setZero();

    Hsparse.resize(nVar, nVar);
    Hsparse.setZero();

    g.resize(nVar);
    g.setZero();

    int rowIdx = 0;
    int colIdx = 0;

    for (const auto &Qi : Q_cost_Matrices)
    {
        // A_cost_all.block(rowIdx, 0, Ai.rows(), Ai.cols()) << Ai;
        EigenUtils::SparseMatrixBlockAsign(Q_cost_all, rowIdx, colIdx, Qi);
        rowIdx += Qi.rows();
        colIdx += Qi.cols();
    }

    rowIdx = 0;
    for (const auto &Ai : A_cost_Matrices)
    {
        // A_cost_all.block(rowIdx, 0, Ai.rows(), Ai.cols()) << Ai;
        EigenUtils::SparseMatrixBlockAsign(A_cost_all, rowIdx, 0, Ai);
        rowIdx += Ai.rows();
    }

    Hsparse += A_cost_all.transpose() * Q_cost_all * A_cost_all;

    rowIdx = 0;
    for (const auto &b_i : b_cost_Vectors)
    {
        b_cost_all.segment(rowIdx, b_i.size()) << b_i;
        rowIdx += b_i.rows();
    }

    g += -A_cost_all.transpose() * Q_cost_all * b_cost_all;
}

// Find the constraint in LinearConstraint_new_regs, add dependency on variables
void MHEproblem::addConstraintDependency(std::string ConstraintName_, std::string VarName, SparseMatrix<double> A)
{
    if (LinearConstraint_new_regs.count(ConstraintName_) && var_idx_regs.count(VarName))
    {
        LinearConstraint_new_regs[ConstraintName_].addDependancy(VarName, A);
    }
    else if (!LinearConstraint_new_regs.count(ConstraintName_))
    {
        std::cout << "Cannot find the constraint " + ConstraintName_ << std::endl;
    }
    else
    {
        std::cout << "Cannot find the variable " + VarName + " for constraint " + ConstraintName_ << std::endl;
    }
}

void MHEproblem::updateConstraintBound(std::string ConstraintName_, VectorXd lb, VectorXd ub, bool equality_update)
{
    if (LinearConstraint_regs.count(ConstraintName_))
    {
        LinearConstraint_regs[ConstraintName_].updateBound(lb, ub, equality_update);
    }
    else
    {
        std::cout << "Cannot find the constraint " + ConstraintName_ << std::endl;
    }
}
void MHEproblem::updateCostGain(std::string CostName_, double scale_gain_)
{
    if (Cost_regs.count(CostName_))
    {
        Cost_regs[CostName_].Q = scale_gain_ * Cost_regs[CostName_].Q;
    }
    else
    {
        std::cout << "Cannot find the Cost " + CostName_ << std::endl;
    }
}
// Find the cost in Cost_new_regs, add dependency on variables
void MHEproblem::addCostDependency(std::string CostName_, std::string VarName, SparseMatrix<double> A)
{
    if (Cost_new_regs.count(CostName_) && var_idx_regs.count(VarName))
    {
        Cost_new_regs[CostName_].addDependancy(VarName, A);
    }
    else if (!Cost_new_regs.count(CostName_))
    {
        std::cout << "cannot find the cost " + CostName_ << std::endl;
    }
    else
    {
        std::cout << "cannot find the variable " + VarName + " for cost " + CostName_ << std::endl;
    }
}

bool MHEproblem::initQP(int T)
{
    // tic("Init");

    bool ifinit = false;

    osqp.clearSolver();

    osqp.data()->setNumberOfConstraints(nConstraints);
    osqp.data()->setNumberOfVariables(nVar);

    osqp.data()->clearHessianMatrix();
    osqp.data()->clearLinearConstraintsMatrix();

    osqp.data()->setHessianMatrix(Hsparse);
    osqp.data()->setGradient(g); // only accepting dense vector
    osqp.data()->setLinearConstraintsMatrix(Aconstrsparse);
    osqp.data()->setLowerBound(lb_all); // only accepting dense vector
    osqp.data()->setUpperBound(ub_all); // only accepting dense vector
    // std::cout << lb_all << std::endl;
    // std::cout << "--------------------------------" << std::endl;
    // std::cout << ub_all << std::endl;

    ifinit = osqp.initSolver();

    // writeToCSVfile("Hsparse.csv", MatrixXd(Hsparse));
    // writeVectorToCSVfile("g.csv", g);
    // writeToCSVfile("Aconstrsparse.csv", MatrixXd(Aconstrsparse));
    // writeVectorToCSVfile("lb_all.csv", lb_all);
    // writeVectorToCSVfile("ub_all.csv", ub_all);

    // update
    // ------------------------------------------------
    // if (T < N + 1)
    // {
    //     osqp.clearSolver();

    //     osqp.data()->setNumberOfConstraints(nConstraints);
    //     osqp.data()->setNumberOfVariables(nVar);

    //     osqp.data()->clearHessianMatrix();
    //     osqp.data()->clearLinearConstraintsMatrix();

    //     osqp.data()->setHessianMatrix(Hsparse);
    //     osqp.data()->setGradient(g);
    //     osqp.data()->setLinearConstraintsMatrix(Aconstrsparse);
    //     osqp.data()->setLowerBound(lb_all);
    //     osqp.data()->setUpperBound(ub_all);

    //     ifinit = osqp.initSolver();
    // }
    // else
    // {
    //     bool up_bound = osqp.updateBounds(lb_all, ub_all);
    //     bool up_gd = osqp.updateGradient(g);

    //     bool up_H = osqp.updateHessianMatrix(Hsparse);
    //     bool up_A = osqp.updateLinearConstraintsMatrix(Aconstrsparse);
    //     // std::cout << up_bound << up_gd << up_H << up_A << std::endl;
    //     ifinit = true;
    // }

    // // Check rank for unconstrained problem
    // JacobiSVD<MatrixXd> svd(MatrixXd(Hsparse));
    // std::cout << "rank: " << svd.rank() << std::endl;

    // SelfAdjointEigenSolver<MatrixXd> solver(H);
    // VectorXd eigen = solver.eigenvalues();
    // std::cout << "Eigen-----------" << std::endl;
    // std::cout << eigen << std::endl;

    // toc("Init");

    return ifinit;
}

void MHEproblem::solveQP()
{

    // tic("Solve");

    osqp.solveProblem();
    solution = osqp.getSolution();

    // toc("Solve");
}

bool MHEproblem::updateQP(int T)
{

    // Todo:
    // 1. add H last block row, last block column.
    // 2. add g last block row.
    // 3. save cost, constraints from new_regs to regs.
    // 4. assemble the cost in Cost_new_regs.
    // 5. add H_new, g_new to H and g.
    // 4. clear cost_new_regs, constraints_new_regs.

    // VO_measurement_string = "VO_measurement_" + std::to_string(T);
    Hsparse.conservativeResize(nVar, nVar);
    EigenUtils::vectorResize(g, nVar);

    // Update QP Cost
    // ------------------------------------
    Hsparse_new.conservativeResize(nVar, nVar);
    Hsparse_new.setZero();

    g_new = VectorXd::Zero(nVar);

    for (const auto &eachcost : Cost_new_regs)
    {
        SparseMatrix<double> Hsparse_cost(nVar, nVar);
        Hsparse_cost.setZero();

        Cost_regs.insert({eachcost.first, eachcost.second});

        SparseMatrix<double> Q_cost = eachcost.second.Q;
        VectorXd b_cost = eachcost.second.b;
        int var_size = b_cost.size();

        for (const auto &cost_Ax_i : eachcost.second.depVarMap)
        {
            std::string var_name_i = cost_Ax_i.first;
            SparseMatrix<double> A_i = cost_Ax_i.second;

            for (const auto &cost_Ax_j : eachcost.second.depVarMap)
            {
                std::string var_name_j = cost_Ax_j.first;
                SparseMatrix<double> A_j = cost_Ax_j.second;
                SparseMatrix<double> H_ij = A_i.transpose() * Q_cost * A_j;

                EigenUtils::SparseMatrixBlockAsign(Hsparse_cost, var_idx_regs[var_name_i] - nVarStart, var_idx_regs[var_name_j] - nVarStart, H_ij);
            }
            VectorXd g_i = -A_i.transpose() * Q_cost * b_cost;

            g_new.segment(var_idx_regs[var_name_i] - nVarStart, g_i.size()) += g_i;
        }
        Hsparse_new += Hsparse_cost;
    }

    Hsparse += Hsparse_new;
    g += g_new;

    Cost_new_regs.clear();
    nCost_new = 0;

    Aconstrsparse.conservativeResize(nConstraints, nVar);
    EigenUtils::vectorResize(lb_all, nConstraints);
    EigenUtils::vectorResize(ub_all, nConstraints);

    Aconstrsparse_new.conservativeResize(nConstraints, nVar);
    Aconstrsparse_new.setZero();

    int row_new = 0;

    // Apending constraints w.r.p to the order they are created
    for (const std::string &constr_name : LinearConstraint_new_name_vec)
    {
        auto eachconstr = LinearConstraint_new_regs[constr_name];

        LinearConstraint_regs.insert({constr_name, eachconstr});

        VectorXd lb_constr = eachconstr.lb;
        VectorXd ub_constr = eachconstr.ub;
        int constr_size = lb_constr.size();
        lb_all.segment(nConstraints - nConstraints_new + row_new, constr_size) << lb_constr;
        ub_all.segment(nConstraints - nConstraints_new + row_new, constr_size) << ub_constr;

        for (const auto &constr_Ax : eachconstr.depVarMap)
        {
            std::string var_name = constr_Ax.first;
            SparseMatrix<double> A_constr = constr_Ax.second;
            EigenUtils::SparseMatrixBlockAsign(Aconstrsparse_new, nConstraints - nConstraints_new + row_new, var_idx_regs[var_name] - nVarStart, A_constr);
        }
        row_new += constr_size;
    }

    Aconstrsparse += Aconstrsparse_new;
    LinearConstraint_new_regs.clear();
    LinearConstraint_new_name_vec.clear();
    nConstraints_new = 0;

    return false;
}

void MHEproblem::Update_Image_bound(std::map<int, Vector3d> vo_constraints_idx_regs, std::map<int, double> vo_reliability_idx_regs)
{
    // VectorXd bound_new = VectorXd::Zero(nConstraints);

    for (const auto &eachvo : vo_constraints_idx_regs)
    {
        lb_all.segment(eachvo.first, 3) << -eachvo.second;
        ub_all.segment(eachvo.first, 3) << -eachvo.second;
        // bound_new.segment(eachvo.first, 3) << eachvo.second;
        // std::cout << "idx" << eachvo.first << std::endl;
    }
    // for (const auto &eachvo : vo_reliability_idx_regs)
    // {
    //     std::string vcam_update_string = "vcam_" + std::to_string(eachvo.first);
    //     // std::cout << vcam_update_string << std::endl;
    //     // std::cout << eachvo.second << std::endl;
    //     // std::cout <<  eachvo.second * Hsparse.block(var_idx_regs[vcam_update_string] - nVarStart, var_idx_regs[vcam_update_string] - nVarStart, dim_cam, dim_cam) << std::endl;
    //     EigenUtils::SparseMatrixBlockAdd(Hsparse, var_idx_regs[vcam_update_string] - nVarStart, var_idx_regs[vcam_update_string] - nVarStart,
    //                                      eachvo.second * Hsparse.block(var_idx_regs[vcam_update_string] - nVarStart, var_idx_regs[vcam_update_string] - nVarStart, dim_cam, dim_cam));
    // }
    // lb_all += bound_new;
    // ub_all += bound_new;
    // std::cout << "n_con" << nConstraints<< std::endl;
    // std::cout << "bound:" << bound_new << std::endl;
}

bool MHEproblem::marginalizeQP(int T)
{
    // get rid of the oldest opt variable
    std::string x_marginalize_string = "x_" + std::to_string(T);
    std::string x_marginalize_leftover_string = "x_" + std::to_string(T + 1);
    std::string v_marginalize_string = "v_" + std::to_string(T);
    std::string w_marginalize_string = "w_" + std::to_string(T);
    std::string vcam_marginalize_string = "vcam_" + std::to_string(T);

    std::string meas_marginalize_string = "Measurement_" + std::to_string(T);
    std::string dyn_marginalize_string = "Dynamic_" + std::to_string(T);
    std::string cam_Meas_marginalize_string = "VO_measurement_" + std::to_string(T);
    std::string contact_marginalize_string = "Contact_" + std::to_string(T);
    std::string lcp_marginalize_string = "Lcp_" + std::to_string(T);
    std::string cost_prior_string = "Prior_0";

    // should have a checking here for safety
    // int var_marginalize_size = dim_state + dim_state + dim_cam + dim_meas;
    // int var_marginalize_size = dim_state + dim_state  + dim_meas;

    int var_marginalize_size = dim_state + dim_state + dim_cam;

    var_idx_regs.erase(x_marginalize_string);
    // var_idx_regs.erase(v_marginalize_string);
    var_idx_regs.erase(w_marginalize_string);
    var_idx_regs.erase(vcam_marginalize_string);

    nVarStart += var_marginalize_size; // update the starting idx of the current window

    nVar = nVarEnd - nVarStart;

    // get rid of the cost correlated to the oldest opt variable
    // Hsparse_new.resize(nVar, nVar);
    // Hsparse_new.setZero();

    Hsparse_new = Hsparse.block(var_marginalize_size, var_marginalize_size,
                                Hsparse.rows() - var_marginalize_size, Hsparse.cols() - var_marginalize_size);

    g_new = g.segment(var_marginalize_size, g.size() - var_marginalize_size);

    if (Cost_regs.count(dyn_marginalize_string))
    {
        // computing the schur complement of x_{i+1}, line by line approach

        // compute the recursive prior M_p and n_p and marginal out H and g
        // 0.5 * (x_0 - \hat(x_0)) ^ T Q_prior_0 (x_0 - \hat(x_0))
        // = 0.5 * (x - b) ^ T Q (x - b) = 0.5 *x ^ T Q x - b ^ T Q x
        if (Cost_regs.count(cost_prior_string))
        {
            M_p = Cost_regs[cost_prior_string].Q;
            n_p = -M_p * Cost_regs[cost_prior_string].b;

            x_prior = Cost_regs[cost_prior_string].b;
            M_p_inv_dense = MatrixXd(M_p).inverse();

            Cost_regs.erase(cost_prior_string);
        }

        solver.compute(M_p);
        SparseMatrixType M_p_inv = solver.solve(I);
        // std::cout << - M_p_inv * n_p<< std::endl;
        if (0)
        {

            SparseMatrix<double> R_meas_marginalize = Cost_regs[meas_marginalize_string].Q;
            SparseMatrix<double> H_meas_marginalize = LinearConstraint_regs[meas_marginalize_string].depVarMap[x_marginalize_string];
            VectorXd y_meas_marginalize = LinearConstraint_regs[meas_marginalize_string].lb;

            SparseMatrix<double> Q_dyn_marginalize = Cost_regs[dyn_marginalize_string].Q;
            SparseMatrix<double> A_dyn_marginalize = LinearConstraint_regs[dyn_marginalize_string].depVarMap[x_marginalize_string];
            VectorXd b_dyn_marginalize = -LinearConstraint_regs[dyn_marginalize_string].lb;

            // solver_dyn.compute(Q_dyn_marginalize);
            // SparseMatrixType Q_dyn_marginalize_inv = solver_dyn.solve(I_dyn);

            // solver_meas.compute(R_meas_marginalize);
            // SparseMatrixType R_meas_marginalize_inv = solver_meas.solve(I_meas);
            MatrixXd Q_inv = MatrixXd(Q_dyn_marginalize).inverse();
            MatrixXd R_inv = MatrixXd(R_meas_marginalize).inverse();
            // std::cout << "----" << std::endl;
            // std::cout << Q_inv << std::endl;
            // std::cout << "----" << std::endl;
            // std::cout << R_inv << std::endl;

            // SparseMatrixType Q_dyn_marginalize_inv = Q_inv.sparseView();
            // SparseMatrixType R_meas_marginalize_inv = R_inv.sparseView();

            // std::cout << M_p_inv_dense << std::endl;
            // std::cout << "+++++++++++++++" << std::endl;

            x_prior = A_dyn_marginalize * x_prior + b_dyn_marginalize;
            M_p_inv_dense = A_dyn_marginalize * M_p_inv_dense * A_dyn_marginalize.transpose() + Q_inv;
            // std::cout << M_p_inv_dense << std::endl;
            // std::cout << "-------" << std::endl;

            MatrixXd temp = MatrixXd(H_meas_marginalize * M_p_inv_dense * H_meas_marginalize.transpose() + R_inv);
            MatrixXd K_KF_ = M_p_inv_dense * H_meas_marginalize.transpose() * temp.inverse();
            x_prior = x_prior + K_KF_ * (y_meas_marginalize - H_meas_marginalize * x_prior);
            M_p_inv_dense = (MatrixXd::Identity(dim_state, dim_state) - K_KF_ * H_meas_marginalize) * M_p_inv_dense;
            M_p = M_p_inv_dense.inverse().sparseView();
            n_p = -M_p * x_prior;
            // std::cout << M_p_inv_dense << std::endl;
            // std::cout << "///////" << std::endl;

            // std::cout << MatrixXd(M_p) << std::endl;
        }
        if (1)
        {
            if (LinearConstraint_regs[cam_Meas_marginalize_string].equality)
            {
                std::cout << "cam" << std::endl;
                SparseMatrix<double> R_meas_marginalize = Cost_regs[meas_marginalize_string].Q;
                SparseMatrix<double> H_meas_marginalize = Cost_regs[meas_marginalize_string].depVarMap[x_marginalize_string];
                VectorXd y_meas_marginalize = Cost_regs[meas_marginalize_string].b;

                SparseMatrix<double> Q_cam_marginalize = Cost_regs[cam_Meas_marginalize_string].Q;
                SparseMatrix<double> A_cam_marginalize = LinearConstraint_regs[cam_Meas_marginalize_string].depVarMap[x_marginalize_string];
                VectorXd b_cam_marginalize = -LinearConstraint_regs[cam_Meas_marginalize_string].lb;

                SparseMatrix<double> Q_dyn_marginalize = Cost_regs[dyn_marginalize_string].Q;
                SparseMatrix<double> A_dyn_marginalize = LinearConstraint_regs[dyn_marginalize_string].depVarMap[x_marginalize_string];
                VectorXd b_dyn_marginalize = -LinearConstraint_regs[dyn_marginalize_string].lb;

                SparseMatrix<double> A_marginalize(dim_cam + dim_state, dim_state);

                EigenUtils::SparseMatrixBlockAsign(A_marginalize, 0, 0, A_dyn_marginalize);
                EigenUtils::SparseMatrixBlockAsign(A_marginalize, dim_state, 0, A_cam_marginalize);

                SparseMatrix<double> Q_marginalize(dim_cam + dim_state, dim_cam + dim_state);

                EigenUtils::SparseMatrixBlockAsign(Q_marginalize, 0, 0, Q_dyn_marginalize);
                EigenUtils::SparseMatrixBlockAsign(Q_marginalize, dim_state, dim_state, Q_cam_marginalize);

                VectorXd b_marginalize = VectorXd::Zero(dim_cam + dim_state);
                b_marginalize.segment(0, dim_state) << b_dyn_marginalize;
                b_marginalize.segment(dim_state, dim_cam) << b_cam_marginalize;

                solver_corelate.compute(Q_marginalize);
                SparseMatrixType Q_marginalize_inv = solver_corelate.solve(I_corelate);

                solver_meas.compute(R_meas_marginalize);
                SparseMatrixType R_meas_marginalize_inv = solver_meas.solve(I_meas);

                SparseMatrix<double> A(dim_state + dim_cam + dim_meas, dim_state + dim_cam + dim_meas);

                A.setZero();

                SparseMatrix<double> A_11 = -A_marginalize * M_p_inv * A_marginalize.transpose();

                A_11 = A_11 - Q_marginalize_inv;

                SparseMatrix<double> A_22 = -H_meas_marginalize * M_p_inv * H_meas_marginalize.transpose();
                A_22 = A_22 - R_meas_marginalize_inv;

                SparseMatrix<double> A_12 = -A_marginalize * M_p_inv * H_meas_marginalize.transpose();
                SparseMatrix<double> A_21 = A_12.transpose();

                EigenUtils::SparseMatrixBlockAsign(A, 0, 0, A_11);
                EigenUtils::SparseMatrixBlockAsign(A, dim_cam + dim_state, dim_cam + dim_state, A_22);
                EigenUtils::SparseMatrixBlockAsign(A, 0, dim_cam + dim_state, A_12);
                EigenUtils::SparseMatrixBlockAsign(A, dim_cam + dim_state, 0, A_21);

                // std::cout << "A" << MatrixXd(A) << std::endl;
                SparseMatrix<double> B(dim_cam + dim_state + dim_meas, dim_state);
                EigenUtils::SparseMatrixBlockAsignFromDense(B, 0, 0, -MatrixXd::Identity(dim_state, dim_state));
                EigenUtils::SparseMatrixBlockAsignFromDense(B, dim_state, 0, -MatrixXd::Identity(dim_cam, dim_cam));

                SparseMatrix<double> C = B.transpose();

                // SparseMatrixType A_inv = MatrixXd(A).llt().solve(MatrixXd::Identity(dim_cam + dim_state + dim_meas,dim_cam + dim_state + dim_meas)).sparseView()); // seems not psotive definte
                SparseMatrixType A_inv = MatrixXd(A).inverse().sparseView();

                VectorXd u = VectorXd::Zero(dim_state + dim_cam + dim_meas);
                u.segment(0, dim_state + dim_cam) = -b_marginalize + A_marginalize * M_p_inv * n_p;
                u.segment(dim_state + dim_cam, dim_meas) = y_meas_marginalize + H_meas_marginalize * M_p_inv * n_p;

                M_p_next = -C * A_inv * B;
                n_p_next = C * A_inv * u;

                M_p = M_p_next;
                n_p = n_p_next;
            }
            else
            {
                SparseMatrix<double> R_meas_marginalize = Cost_regs[meas_marginalize_string].Q;
                SparseMatrix<double> H_meas_marginalize = Cost_regs[meas_marginalize_string].depVarMap[x_marginalize_string];
                VectorXd y_meas_marginalize = Cost_regs[meas_marginalize_string].b;

                SparseMatrix<double> Q_dyn_marginalize = Cost_regs[dyn_marginalize_string].Q;
                SparseMatrix<double> A_dyn_marginalize = LinearConstraint_regs[dyn_marginalize_string].depVarMap[x_marginalize_string];
                VectorXd b_dyn_marginalize = -LinearConstraint_regs[dyn_marginalize_string].lb;

                solver_dyn.compute(Q_dyn_marginalize);
                SparseMatrixType Q_dyn_marginalize_inv = solver_dyn.solve(I_dyn);

                solver_meas.compute(R_meas_marginalize);
                SparseMatrixType R_meas_marginalize_inv = solver_meas.solve(I_meas);

                SparseMatrix<double> A(dim_state + dim_meas, dim_state + dim_meas);

                A.setZero();

                SparseMatrix<double> A_11 = -A_dyn_marginalize * M_p_inv * A_dyn_marginalize.transpose();

                A_11 = A_11 - Q_dyn_marginalize_inv;

                SparseMatrix<double> A_22 = -H_meas_marginalize * M_p_inv * H_meas_marginalize.transpose();
                A_22 = A_22 - R_meas_marginalize_inv;

                SparseMatrix<double> A_12 = -A_dyn_marginalize * M_p_inv * H_meas_marginalize.transpose();
                SparseMatrix<double> A_21 = A_12.transpose();

                EigenUtils::SparseMatrixBlockAsign(A, 0, 0, A_11);
                EigenUtils::SparseMatrixBlockAsign(A, dim_state, dim_state, A_22);
                EigenUtils::SparseMatrixBlockAsign(A, 0, dim_state, A_12);
                EigenUtils::SparseMatrixBlockAsign(A, dim_state, 0, A_21);

                SparseMatrix<double> B(dim_state + dim_meas, dim_state);
                EigenUtils::SparseMatrixBlockAsignFromDense(B, 0, 0, -MatrixXd::Identity(dim_state, dim_state));

                SparseMatrix<double> C = B.transpose();

                SparseMatrixType A_inv = MatrixXd(A).inverse().sparseView();

                VectorXd u = VectorXd::Zero(dim_state + dim_meas);
                u.segment(0, dim_state) = -b_dyn_marginalize + A_dyn_marginalize * M_p_inv * n_p;
                u.segment(dim_state, dim_meas) = y_meas_marginalize + H_meas_marginalize * M_p_inv * n_p;

                M_p_next = -C * A_inv * B;
                n_p_next = C * A_inv * u;

                M_p = M_p_next;
                n_p = n_p_next;

                // std::cout << "prior" << std::endl;

                // std::cout << (-MatrixXd(M_p).inverse() * n_p).segment<12>(27) << std::endl;
            }
        }
        // if (1)
        // {
        //     SparseMatrix<double> R_meas_marginalize = Cost_regs[meas_marginalize_string].Q;
        //     SparseMatrix<double> H_meas_marginalize = LinearConstraint_regs[meas_marginalize_string].depVarMap[x_marginalize_string];
        //     VectorXd y_meas_marginalize = LinearConstraint_regs[meas_marginalize_string].lb;

        //     // SparseMatrix<double> R_meas_marginalize = Cost_regs[meas_marginalize_string].Q;
        //     // SparseMatrix<double> H_meas_marginalize = Cost_regs[meas_marginalize_string].depVarMap[x_marginalize_string];
        //     // VectorXd y_meas_marginalize = Cost_regs[meas_marginalize_string].b;

        //     SparseMatrix<double> Q_dyn_marginalize = Cost_regs[dyn_marginalize_string].Q;
        //     SparseMatrix<double> A_dyn_marginalize = LinearConstraint_regs[dyn_marginalize_string].depVarMap[x_marginalize_string];
        //     VectorXd b_dyn_marginalize = LinearConstraint_regs[dyn_marginalize_string].lb;

        //     SparseMatrix<double> A_constraint_marginalize = LinearConstraint_regs[contact_marginalize_string].depVarMap[x_marginalize_string];
        //     VectorXd b_constraint_marginalize = LinearConstraint_regs[contact_marginalize_string].lb;
        //     int dim_contact = b_constraint_marginalize.size();

        //     solver_dyn.compute(Q_dyn_marginalize);
        //     SparseMatrixType Q_dyn_marginalize_inv = solver_dyn.solve(I_dyn);
        //     solver_meas.compute(R_meas_marginalize);
        //     SparseMatrixType R_meas_marginalize_inv = solver_meas.solve(I_meas);

        //     SparseMatrix<double> A(2 * dim_state + dim_meas, 2 * dim_state + dim_meas);
        //     A.setZero();
        //     EigenUtils::SparseMatrixBlockAsign(A, 0, 0, M_p_inv);
        //     EigenUtils::SparseMatrixBlockAsign(A, dim_state, dim_state, R_meas_marginalize_inv);
        //     EigenUtils::SparseMatrixBlockAsign(A, dim_state + dim_meas, dim_state + dim_meas, Q_dyn_marginalize_inv);

        //     // SparseMatrix<double> C(dim_meas + dim_state + dim_contact, 2 * dim_state + dim_meas);
        //     SparseMatrix<double> C(dim_meas + dim_state, 2 * dim_state + dim_meas);
        //     C.setZero();
        //     EigenUtils::SparseMatrixBlockAsign(C, 0, 0, H_meas_marginalize);
        //     EigenUtils::SparseMatrixBlockAsignFromDense(C, 0, dim_state, MatrixXd::Identity(dim_meas, dim_meas));
        //     EigenUtils::SparseMatrixBlockAsign(C, dim_meas, 0, A_dyn_marginalize);
        //     EigenUtils::SparseMatrixBlockAsignFromDense(C, dim_meas, dim_state + dim_meas, MatrixXd::Identity(dim_state, dim_state));
        //     // EigenUtils::SparseMatrixBlockAsign(C, dim_meas + dim_state, 0, A_constraint_marginalize);

        //     SparseMatrix<double> D = -C * A * C.transpose();

        //     // std::cout << MatrixXd(A_constraint_marginalize) << std::endl;

        //     // SparseMatrix<double> D_12(dim_state + dim_meas + dim_contact, dim_state);
        //     SparseMatrix<double> D_12(dim_state + dim_meas, dim_state);
        //     EigenUtils::SparseMatrixBlockAsignFromDense(D_12, dim_meas, 0, -MatrixXd::Identity(dim_state, dim_state));

        //     SparseMatrix<double> D_21 = D_12.transpose();

        //     SparseMatrixType D_inv = MatrixXd(D).inverse().sparseView();
        //     // std::cout << MatrixXd(D_inv) << std::endl;

        //     VectorXd u_up = VectorXd::Zero(2 * dim_state + dim_meas);
        //     u_up.segment(0, dim_state) = -n_p;

        //     VectorXd u_down = VectorXd::Zero(dim_state + dim_meas + dim_contact);
        //     // VectorXd u_down = VectorXd::Zero(dim_state + dim_meas);
        //     u_down.segment(0, dim_meas) = y_meas_marginalize;
        //     u_down.segment(dim_meas, dim_state) = b_dyn_marginalize;
        //     u_down.segment(dim_meas + dim_state, dim_contact) = b_constraint_marginalize;

        //     M_p_next = -D_21 * D_inv * D_12;
        //     n_p_next = D_21 * D_inv * (u_down - C * A * u_up);

        //     M_p = M_p_next;
        //     n_p = n_p_next;

        //     // // ------------------------------------
        //     // // check distribution
        //     // std::cout << "Diagonal elements: ";
        //     // MatrixXd M_p_dense = MatrixXd(M_p);
        //     // for (int i = 0; i < M_p_dense.rows(); ++i)
        //     // {
        //     //     std::cout << M_p_dense(i, i) << " ";
        //     // }
        //     // std::cout << std::endl;
        // }

        // combined terms
        SparseMatrix<double> H_marginalize_sparse = M_p;
        VectorXd g_marginalize = n_p;

        // Resize to nVar, nVar to operate with Hsparse
        H_marginalize_sparse.conservativeResize(nVar, nVar);
        // g_marginalize.conservativeResize(nVar);
        EigenUtils::vectorResize(g_marginalize, nVar);

        Hsparse.resize(nVar, nVar);
        Hsparse.setZero();
        g.resize(nVar);
        g.setZero();

        Hsparse += Hsparse_new + H_marginalize_sparse;
        g += g_new + g_marginalize;

        nCost += -Cost_regs[dyn_marginalize_string].b.size();
        nCost += -Cost_regs[cam_Meas_marginalize_string].b.size();
        nCost += -Cost_regs[meas_marginalize_string].b.size();

        Cost_regs.erase(dyn_marginalize_string);
        Cost_regs.erase(cam_Meas_marginalize_string);
        Cost_regs.erase(meas_marginalize_string);
    }
    else
    {
        std::cout << "Error missing dynamic cost at:" + std::to_string(T) << std::endl;
    }
    if (!(nConstraints == 0))
    {
        int constraints_marginalize_size = 0;
        constraints_marginalize_size += LinearConstraint_regs[dyn_marginalize_string].lb.size();
        constraints_marginalize_size += LinearConstraint_regs[cam_Meas_marginalize_string].lb.size();
        // constraints_marginalize_size += LinearConstraint_regs[meas_marginalize_string].lb.size();
        constraints_marginalize_size += LinearConstraint_regs[contact_marginalize_string].lb.size();
        // constraints_marginalize_size += LinearConstraint_regs[lcp_marginalize_string].lb.size();

        nConstraints += -constraints_marginalize_size;

        // Aconstrsparse_new.resize(nConstraints, nVar);
        Aconstrsparse_new = Aconstrsparse.block(constraints_marginalize_size, var_marginalize_size, Aconstrsparse.rows() - constraints_marginalize_size, Aconstrsparse.cols() - var_marginalize_size);
        lb_all_new = lb_all.segment(constraints_marginalize_size, lb_all.size() - constraints_marginalize_size);
        ub_all_new = ub_all.segment(constraints_marginalize_size, ub_all.size() - constraints_marginalize_size);

        Aconstrsparse.resize(nConstraints, nVar);
        Aconstrsparse.setZero();
        Aconstrsparse += Aconstrsparse_new;

        lb_all.resize(nConstraints);
        lb_all.setZero();
        lb_all += lb_all_new;

        ub_all.resize(nConstraints);
        ub_all.setZero();
        ub_all += ub_all_new;

        LinearConstraint_regs.erase(dyn_marginalize_string);
        LinearConstraint_regs.erase(cam_Meas_marginalize_string);
        // LinearConstraint_regs.erase(meas_marginalize_string);
        LinearConstraint_regs.erase(contact_marginalize_string);
        // LinearConstraint_regs.erase(lcp_marginalize_string);
    }
    return false;
}

void MHEproblem::getsolution(int T)
{
    std::string x_name = "x_" + std::to_string(T);
    std::string w_name = "w_" + std::to_string(T - 1);
    std::string v_name = "v_" + std::to_string(T);
    std::string vcam_name = "vcam_" + std::to_string(T - 1);
    x_solution_now = solution.segment(var_idx_regs[x_name] - nVarStart, dim_state);
    // std::cout << "v_res" << std::endl;
    // std::cout << solution.segment(var_idx_regs[v_name] - nVarStart, dim_meas) << std::endl;
    // std::cout << "w_res" << std::endl;
    // std::cout << solution.segment(var_idx_regs[w_name] - nVarStart, dim_state) << std::endl;
}

// Formulate optimization problem at once from Cost_regs and Constraint_regs
bool MHEproblem::formulateQP()
{
    // formulate and initialize the QP and returns if success
    assembleCost();
    assembleConstraints();
    return true;
}

void MHEproblem::resetQP()
{
    // clear the solver
    osqp.clearSolver();

    // clear all elements
    nVar = 0;
    nConstraints = 0;
    nCost = 0;

    // opt varable
    var_idx_regs.clear();

    // linear constraints
    LinearConstraint_regs.clear();

    A_constraint_Matrices.clear();
    lb_vectors.clear();
    ub_vectors.clear();

    // cost: norm of linear eqs.
    Cost_regs.clear();
    Cost_new_regs.clear();
    A_cost_Matrices.clear();
    b_cost_Vectors.clear();
    Q_cost_Matrices.clear();
}

void MHEproblem::tic(std::string str, int mode)
{
    static std::chrono::_V2::system_clock::time_point t_start;

    if (mode == 0)
        t_start = std::chrono::high_resolution_clock::now();
    else
    {
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout << str + " elapsed time: " << (t_end - t_start).count() * 1E-9 << " seconds\n";
    }
}

void MHEproblem::toc(std::string str) { tic(str, 1); }
