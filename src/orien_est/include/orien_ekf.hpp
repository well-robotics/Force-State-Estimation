// IMU Quaternion EKF, implemented based on https://ahrs.readthedocs.io/en/latest/filters/ekf.html

#ifndef IMU_EKF_HPP
#define IMU_EKF_HPP

#include <rclcpp/rclcpp.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <sensor_msgs/msg/imu.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <optitrack_broadcast/msg/mocap.hpp>
#include "std_msgs/msg/header.hpp"

// ros msg
#include "unitree_go/msg/low_state.hpp"
#include "unitree_go/msg/sport_mode_state.hpp"

using namespace Eigen;

namespace orien_ekf
{

    class orien_ekf : public rclcpp::Node
    {
    private:
        rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub;
        void vo_pose_callback(const geometry_msgs::msg::PoseStamped::SharedPtr msg);

        rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr vo_pose_sub;
        void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg);

        // Ground Truth
        rclcpp::Subscription<optitrack_broadcast::msg::Mocap>::SharedPtr mocap_sub;
        void mocap_callback(const optitrack_broadcast::msg::Mocap::SharedPtr msg);

        rclcpp::Subscription<unitree_go::msg::LowState>::SharedPtr low_state_sub;
        void low_state_callback(const unitree_go::msg::LowState::SharedPtr msg);

        rclcpp::TimerBase::SharedPtr timer_;
        void timerCallback();

        rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr publisher_filter_;

        double dt_ = 0.002;

        Vector3d gravity_ = Vector3d::Zero();

        double time_init_;
        double imu_time_;
        double vo_time_;
        int discrete_time_ = 0;
        std_msgs::msg::Header imu_header_;

        int using_vo;

    private:
        Vector3d accel_b_ = Vector3d::Zero();
        Vector3d angular_vel_b_ = Vector3d::Zero();

        bool vo_new_ = false;
        VectorXd vo_pose_quaternion_ = VectorXd::Zero(4);
        VectorXd vo_off_ = VectorXd::Zero(4);

        std::vector<Vector3d> angular_vel_stack_;
        std::vector<Vector3d> accel_stack_;
        std::vector<double> imu_time_stack_;
        std::vector<int> discrete_time_stack_;

        std::vector<Vector4d> filter_quaternion_stack_;
        std::vector<Matrix4d> filter_Cov_stack_;

    private:
        VectorXd quaternion_ = VectorXd::Zero(4);
        VectorXd gt_quaternion_ = VectorXd::Zero(4);

        VectorXd quaternion_pred_ = VectorXd::Zero(4);
        VectorXd quaternion_correct_ = VectorXd::Zero(4);
        MatrixXd Cov_q_ = MatrixXd::Zero(4, 4);
        MatrixXd Cov_q_pred_ = MatrixXd::Zero(4, 4);
        MatrixXd Cov_q_correct_ = MatrixXd::Zero(4, 4);

        Matrix3d C_accel_ = Matrix3d::Zero();
        Matrix3d C_gyro_ = Matrix3d::Zero();
        Matrix4d C_vo_ = MatrixXd::Zero(4, 4);

        int init_vo = 0;
        int init_imu = 0;
        int init_filter = 0;
        int init_cali = 0;

    public:
        orien_ekf(const std::string &name);

        void gyro_nonlinear_predict(VectorXd &quaternion_pred, VectorXd &quaternion, Vector3d &gyro_reading, MatrixXd &Cov_q, MatrixXd &Cov_q_pred);
        void gyro_nonlinear_correct(VectorXd &quaternion_correct, VectorXd &quaternion_pred, Vector3d &accel_readings, MatrixXd &Cov_q_pred, MatrixXd &Cov_q_correct);
        void vo_nonlinear_correct(VectorXd &quaternion_correct, VectorXd &quaternion_pred, VectorXd &quaternion_vo, MatrixXd &Cov_q_pred, MatrixXd &Cov_q_correct);

        void get_measurement();

        void gyro_2_Ohm(Vector3d &gyro_reading, MatrixXd &Ohm);
        void quat_2_W(VectorXd &quaternion, MatrixXd &W);
        void quat_2_Rot(VectorXd &quaternion, MatrixXd &Rotation);
        void quat_2_H(VectorXd &quaternion, MatrixXd &H);
        void quat_norm(VectorXd &quaternion);
        void quaternionToEuler(VectorXd quaternion, double &yaw, double &pitch, double &roll);
        void quat_mul(VectorXd &quaternion_left, VectorXd &quaternion_right, VectorXd &quaternion_output);
        void quat_inv(VectorXd &quaternion, VectorXd &quaternion_out);
    };
}
#endif // NGIMU_TF_HPP