#ifndef EIGEN_UTILS_H
#define EIGEN_UTILS_H
#include <Eigen/Sparse>

namespace EigenUtils
{
    using namespace Eigen;

    inline void SparseMatrixBlockAsign(SparseMatrix<double> &Aout, int i, int j, SparseMatrix<double> Ain)
    {
        // block assignment of sparse matrix. (didn't find something like this in Eigen, implement it here temporarly)
        int idxRow, idxCol;
        for (int k = 0; k < Ain.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(Ain, k); it; ++it)
            {
                idxRow = i + it.row();
                idxCol = j + it.col();
                Aout.insert(idxRow, idxCol) = it.value();
                // Aout.coeffRef(idxRow, idxCol) = it.value();
            }
    }

    inline void SparseMatrixBlockAdd(SparseMatrix<double> &Aout, int i, int j, SparseMatrix<double> Ain)
    {
        // block assignment of sparse matrix. (didn't find something like this in Eigen, implement it here temporarly)
        int idxRow, idxCol;
        for (int k = 0; k < Ain.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(Ain, k); it; ++it)
            {
                idxRow = i + it.row();
                idxCol = j + it.col();
                Aout.coeffRef(idxRow, idxCol) += it.value();
            }
    }

    inline void SparseMatrixBlockReplaceFromDense(SparseMatrix<double> &Aout, int i, int j, MatrixXd Ain)
    {
        // block assignment of sparse matrix. (didn't find something like this in Eigen, implement it here temporarly)
        int idxRow, idxCol;
        for (int row = 0; row < Ain.rows(); ++row)
        {
            for (int col = 0; col < Ain.cols(); ++col)
            {
                idxRow = i + row;
                idxCol = j + col;
                Aout.coeffRef(idxRow, idxCol) = Ain(row, col);
            }
        }
    }

    inline void SparseMatrixBlockAsignFromDense(SparseMatrix<double> &Aout, int i, int j, MatrixXd Ain)
    {
        // block assignment of sparse matrix. (didn't find something like this in Eigen, implement it here temporarly)
        int idxRow, idxCol;
        for (int row = 0; row < Ain.rows(); ++row)
        {
            for (int col = 0; col < Ain.cols(); ++col)
            {
                idxRow = i + row;
                idxCol = j + col;
                Aout.insert(idxRow, idxCol) = Ain(row, col);
            }
        }
    }

    inline void SparseMatrixBlockAcquisition(SparseMatrix<double> &Aout, int i, int j, SparseMatrix<double> Ain)
    {
        // // block acquisition of sparse matrix. (almost the same behavior compared with .block() operation)
        int idxRow = 0;
        int idxCol = 0;
        for (int k = j; k < Ain.outerSize(); ++k)
        {
            for (SparseMatrix<double>::InnerIterator it(Ain, k); it; ++it)
            {

                if (it.row() >= i)
                {
                    idxRow = it.row() - i;
                    idxCol = it.col() - j;
                }
                Aout.insert(idxRow, idxCol) = it.value();
                // Aout.coeffRef(idxRow, idxCol) = it.value();
            }
        }
        // std::vector<Triplet<double>> tripletList;
        // for (int k = i; k < Ain.outerSize(); ++k)
        // {
        //     for (SparseMatrix<double>::InnerIterator it(Ain, k); it; ++it)
        //     {
        //         if (it.row() >= j)
        //         {
        //             tripletList.push_back(Triplet<double>(it.row() - i, it.col() - j, it.value()));
        //         }
        //     }
        // }
        // Aout.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    inline void vectorResize(VectorXd &vector, const int new_size)
    {
        int old_size = vector.size();
        vector.conservativeResize(new_size);
        vector.segment(old_size, new_size - old_size).setZero();
    }

    inline void vector3dSkew(Matrix3d &skew_sym, const Vector3d &vector)
    {

        skew_sym << 0, -vector(2), vector(1),
            vector(2), 0, -vector(0),
            -vector(1), vector(0), 0;
    }

    inline void QuaternionToEuler(Quaterniond quaternion, Vector3d &euler)
    {
        double roll, pitch, yaw = 0.0;
        double qw = quaternion.w();
        double qx = quaternion.x();
        double qy = quaternion.y();
        double qz = quaternion.z();
        // roll (x-axis rotation)
        double sinr_cosp = 2 * (qw * qx + qy * qz);
        double cosr_cosp = 1 - 2 * (qx * qx + qy * qy);
        roll = atan2(sinr_cosp, cosr_cosp);
        // pitch (y-axis rotation)
        double sinp = 2 * (qw * qy - qz * qx);
        if (abs(sinp) >= 1)
            pitch = copysign(M_PI / 2, sinp); // Use M_PI for pi in C++
        else
            pitch = asin(sinp);
        // yaw (z-axis rotation)
        double siny_cosp = 2 * (qw * qz + qx * qy);
        double cosy_cosp = 1 - 2 * (qy * qy + qz * qz);
        yaw = atan2(siny_cosp, cosy_cosp);
        euler(0) = roll;
        euler(1) = pitch;
        euler(2) = yaw;
    }

    inline void VectorToMatrixXd(const std::vector<std::vector<double>> &vec, MatrixXd &mat)
    {
        if (vec.empty() || vec[0].empty())
        {
            throw std::invalid_argument("Input vector is empty or improperly formed.");
        }

        size_t rows = vec.size();
        size_t cols = vec[0].size();

        mat = MatrixXd::Zero(rows, cols);

        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < cols; ++j)
            {
                mat(i, j) = vec[i][j];
            }
        }
    }

    inline void AngularToEulerRates(const Vector3d &angular_velocity, const Vector3d &euler_angles, Vector3d &euler_rates)
    {
        double phi = euler_angles(0);   // Roll (phi)
        double theta = euler_angles(1); // Pitch (theta)
        double psi = euler_angles(2);   // Yaw (psi)

        // Create the transformation matrix from angular velocity to Euler rates
        Matrix3d T;
        T << 1, std::sin(phi) * std::tan(theta), std::cos(phi) * std::tan(theta),
            0, std::cos(phi), -std::sin(phi),
            0, std::sin(phi) / std::cos(theta), std::cos(phi) / std::cos(theta);

        // Calculate the Euler rates
        euler_rates = T * angular_velocity;
    }

    inline void deltaEuler2deltaQuaternion(const Vector3d &deltaEuler, Quaterniond &deltaQuat)
    {
        // Extract individual delta angles
        double dRoll = deltaEuler(0);  // Rotation around X axis
        double dPitch = deltaEuler(1); // Rotation around Y axis
        double dYaw = deltaEuler(2);   // Rotation around Z axis

        // Delta quaternion for roll (X axis)
        Eigen::Quaterniond dqX(std::cos(dRoll / 2), std::sin(dRoll / 2), 0.0, 0.0);

        // Delta quaternion for pitch (Y axis)
        Eigen::Quaterniond dqY(std::cos(dPitch / 2), 0.0, std::sin(dPitch / 2), 0.0);

        // Delta quaternion for yaw (Z axis)
        Eigen::Quaterniond dqZ(std::cos(dYaw / 2), 0.0, 0.0, std::sin(dYaw / 2));

        // Combine the quaternions in ZYX order (Z * Y * X)
        deltaQuat = dqZ * dqY * dqX;
    }

}
#endif // EIGEN_UTILS_H