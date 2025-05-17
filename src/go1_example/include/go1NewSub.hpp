#ifndef GO1_NEW_SUB_HPP
#define GO1_NEW_SUB_HPP

#include <rclcpp/rclcpp.hpp>

// CPP header(non ros)
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include "decentral_legged_est/EigenUtils.hpp"
#include "decentral_legged_est/EstSub.hpp"

// ros msg
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/joint_state.hpp>
#include <std_msgs/msg/int16_multi_array.hpp>
#include <optitrack_broadcast/msg/mocap.hpp>
#include <nav_msgs/msg/odometry.hpp>

// kinematics lib
#include "../src/Expressions/FL_foot.hh"
#include "../src/Expressions/FR_foot.hh"
#include "../src/Expressions/RL_foot.hh"
#include "../src/Expressions/RR_foot.hh"
#include "../src/Expressions/J_FL.hh"
#include "../src/Expressions/J_FR.hh"
#include "../src/Expressions/J_RL.hh"
#include "../src/Expressions/J_RR.hh"

using namespace Eigen;

namespace robotSub
{

    class go1NewSub : public robotSub
    {
    private:
        // IMU
        rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub;
        void imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg);

        // Leg Odometry, encoder information & kinematics
        rclcpp::Subscription<sensor_msgs::msg::JointState>::SharedPtr lo_sub;
        void lo_callback(const sensor_msgs::msg::JointState::SharedPtr msg);

        // Contact
        rclcpp::Subscription<std_msgs::msg::Int16MultiArray>::SharedPtr contact_sub;
        void contact_callback(const std_msgs::msg::Int16MultiArray::SharedPtr msg);

        // Ground Truth
        rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr mocap_sub;
        void mocap_callback(const nav_msgs::msg::Odometry::SharedPtr msg);

    public:
        go1NewSub(const std::string &name);
        ~go1NewSub();
    };
}
#endif // GO1_NEW_SUB_HPP