#include "go1NewSub.hpp"

using namespace Eigen;

namespace robotSub
{

    go1NewSub::go1NewSub(const std::string &name) : robotSub(name)
    {

        // Go1 sensors
        // IMU
        imu_sub = create_subscription<sensor_msgs::msg::Imu>("/b1/imu",
                                                             10,
                                                             std::bind(&go1NewSub::imu_callback, this, std::placeholders::_1));
        // Leg Odometry, encoder information & kinematics
        lo_sub = create_subscription<sensor_msgs::msg::JointState>("/b1/joint_states",
                                                                   10,
                                                                   std::bind(&go1NewSub::lo_callback, this, std::placeholders::_1));

        // Contact
        contact_sub = create_subscription<std_msgs::msg::Int16MultiArray>("/b1/contact",
                                                                          10,
                                                                          std::bind(&go1NewSub::contact_callback, this, std::placeholders::_1));
        // Ground Truth
        mocap_sub = create_subscription<nav_msgs::msg::Odometry>("/b1/onboard",
                                                                 10,
                                                                 std::bind(&go1NewSub::mocap_callback, this, std::placeholders::_1));

    }

    go1NewSub::~go1NewSub()
    {
    }

    void go1NewSub::contact_callback(const std_msgs::msg::Int16MultiArray::SharedPtr msg)
    {
        //---------------------------------------------------------------
        // ToDo: configurate/fill-up the contact;
        // robot_store_->contact_: VectorXd (num_legs) (0: no_contact, 1: contact);
        //---------------------------------------------------------------
        VectorXd contact = VectorXd::Zero(robot_params_->num_legs_);

        for (int i = 0; i < robot_params_->num_legs_; i++)
        {
            contact(i) = (static_cast<double>(msg->data[i]) >= robot_params_->contact_effort_theshold_) ? 1.0 : 0.0;
        }
        robot_store_->contact_ = contact;
    }

    void go1NewSub::imu_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
    {
        //---------------------------------------------------------------
        // ToDo: configurate/fill-up the imu_time_, accel_b_, angular_b_;
        // robot_store_->imu_time_: double;
        // robot_store_->accel_b_: Vector3d;
        // robot_store_->angular_b_: Vector3d;
        //---------------------------------------------------------------

        rclcpp::Time imu_time = rclcpp::Clock().now();
        robot_store_->imu_time_ = static_cast<double>(imu_time.nanoseconds()) / 1e9 - time_init_; // Correctly obtaining the timestamp in seconds as a double


        robot_store_->accel_b_(0) = msg->linear_acceleration.x;
        robot_store_->accel_b_(1) = msg->linear_acceleration.y;
        robot_store_->accel_b_(2) = msg->linear_acceleration.z;

        robot_store_->angular_b_(0) = msg->angular_velocity.x;
        robot_store_->angular_b_(1) = msg->angular_velocity.y;
        robot_store_->angular_b_(2) = msg->angular_velocity.z;

        Quaterniond imu_quat;
        imu_quat.x() = msg->orientation.x;
        imu_quat.y() = msg->orientation.y;
        imu_quat.z() = msg->orientation.z;
        imu_quat.w() = msg->orientation.w;

        robot_store_->quaternion_ = imu_quat;

        // // Test with Unitree IMU orientation
        // robot_store_->quaternion_.x() = msg->orientation.x;
        // robot_store_->quaternion_.y() = msg->orientation.y;
        // robot_store_->quaternion_.z() = msg->orientation.z;
        // robot_store_->quaternion_.w() = msg->orientation.w;

        // Quaterniond quat_off = imu_quat * delta_quat;
        // std::cout << delta_euler << std::endl;
        // std::cout << delta_quat.x() << "," << delta_quat.y() << "," << delta_quat.z() << "," << delta_quat.w() << std::endl;

        // std::cout << quat_off.x() << "," << quat_off.y() << "," << quat_off.z() << "," << quat_off.w() << std::endl;

        EigenUtils::QuaternionToEuler(robot_store_->quaternion_, filter_euler_);

        imu_msg_num_++;
    }

    void go1NewSub::lo_callback(const sensor_msgs::msg::JointState::SharedPtr msg)
    {
        //---------------------------------------------------------------
        // ToDo: configurate/fill-up the encoder_position, encoder_velocity, foot_forward_kinematics, foot_jacobian, contact;
        // robot_store_->joint_states_position_: VectorXd(num_joints*num_legs);
        // robot_store_->joint_states_velocity_: VectorXd(num_joints*num_legs);
        // robot_store_->p_imu_2_foot_: MatrixXd, (3*num_legs,1);
        // robot_store_->J_imu_2_foot_: MatrixXd, (3*num_legs,3);
        //---------------------------------------------------------------
        VectorXd joint_position = Map<VectorXd, Unaligned>(const_cast<double *>(msg->position.data()), msg->position.size());
        VectorXd joint_velocity = Map<VectorXd, Unaligned>(const_cast<double *>(msg->velocity.data()), msg->velocity.size());
        VectorXd joint_effort = Map<VectorXd, Unaligned>(const_cast<double *>(msg->effort.data()), msg->effort.size());

        robot_store_->joint_states_position_ = joint_position.segment(0, 3 * robot_params_->num_legs_);
        robot_store_->joint_states_velocity_ = joint_velocity.segment(0, 3 * robot_params_->num_legs_);
        robot_store_->joint_states_effort_ = joint_effort.segment(0, 3 * robot_params_->num_legs_);

        VectorXd joint_position_append = VectorXd::Zero(18); // [p,v, 4 foot * 3 joints, ]
        VectorXd joint_velocity_append = VectorXd::Zero(18); // [p,v, 4 foot * 3 joints, ]

        joint_position_append.segment(6, 3 * robot_params_->num_legs_) = joint_position.segment(0, 3 * robot_params_->num_legs_);
        joint_velocity_append.segment(6, 3 * robot_params_->num_legs_) = joint_velocity.segment(0, 3 * robot_params_->num_legs_);

        MatrixXd p_imu_2_foot_ = MatrixXd::Zero(12, 1);
        MatrixXd J_imu_2_foot_ = MatrixXd::Zero(12, 3);

        MatrixXd p_ib = MatrixXd::Zero(1, 3);
        p_ib << robot_params_->p_ib_[0], robot_params_->p_ib_[1], robot_params_->p_ib_[2]; // Modeled in body frame (IMU twisted but not moved, so p_ib_b is fixed)

        Quaterniond q_ib;
        q_ib.w() = robot_params_->quaternion_ib_[0];
        q_ib.x() = robot_params_->quaternion_ib_[1];
        q_ib.y() = robot_params_->quaternion_ib_[2];
        q_ib.z() = robot_params_->quaternion_ib_[3];

        // Matrix3d R_ib = q_ib.normalized().toRotationMatrix();
        Matrix3d R_ib = Matrix3d::Identity();
        // Foot order: FR, FL, RR, RL; for both hardware and kinematics lib
        // Joints Order: hip, thig, calf, foot(fixed 0); for both hardware (no foot joint) and kinematics lib
        MatrixXd p_FR_foot = MatrixXd::Zero(1, 3);
        SymFunction::FR_foot(p_FR_foot, joint_position_append);
        p_FR_foot = (p_FR_foot + p_ib);
        p_imu_2_foot_.block<3, 1>(0 * 3, 0) = p_FR_foot.transpose();

        MatrixXd J_FR_foot = MatrixXd::Zero(3, 18);
        SymFunction::J_FR(J_FR_foot, joint_position_append);
        J_imu_2_foot_.block<3, 3>(0 * 3, 0) = J_FR_foot.block<3, 3>(0, 6 + 0 * 3);

        // FL
        MatrixXd p_FL_foot = MatrixXd::Zero(1, 3);
        SymFunction::FL_foot(p_FL_foot, joint_position_append);
        p_FL_foot = (p_FL_foot + p_ib);
        p_imu_2_foot_.block<3, 1>(1 * 3, 0) = p_FL_foot.transpose();

        MatrixXd J_FL_foot = MatrixXd::Zero(3, 18);
        SymFunction::J_FL(J_FL_foot, joint_position_append);
        J_imu_2_foot_.block<3, 3>(1 * 3, 0) = J_FL_foot.block<3, 3>(0, 6 + 1 * 3);

        // RR
        MatrixXd p_RR_foot = MatrixXd::Zero(1, 3);
        SymFunction::RR_foot(p_RR_foot, joint_position_append);
        p_RR_foot = (p_RR_foot + p_ib);
        p_imu_2_foot_.block<3, 1>(2 * 3, 0) = p_RR_foot.transpose();

        MatrixXd J_RR_foot = MatrixXd::Zero(3, 18);
        SymFunction::J_RR(J_RR_foot, joint_position_append);
        J_imu_2_foot_.block<3, 3>(2 * 3, 0) = J_RR_foot.block<3, 3>(0, 6 + 2 * 3);

        // RL
        MatrixXd p_RL_foot = MatrixXd::Zero(1, 3);
        SymFunction::RL_foot(p_RL_foot, joint_position_append);
        p_RL_foot = (p_RL_foot + p_ib);
        p_imu_2_foot_.block<3, 1>(3 * 3, 0) = p_RL_foot.transpose();

        MatrixXd J_RL_foot = MatrixXd::Zero(3, 18);
        SymFunction::J_RL(J_RL_foot, joint_position_append);
        J_imu_2_foot_.block<3, 3>(3 * 3, 0) = J_RL_foot.block<3, 3>(0, 6 + 3 * 3);

        robot_store_->p_imu_2_foot_ = p_imu_2_foot_;
        robot_store_->J_imu_2_foot_ = J_imu_2_foot_;
    }

    void go1NewSub::mocap_callback(const nav_msgs::msg::Odometry::SharedPtr msg)
    {

        robot_store_->gt_p_(0) = msg->pose.pose.position.x;
        robot_store_->gt_p_(1) = msg->pose.pose.position.y;
        robot_store_->gt_p_(2) = msg->pose.pose.position.z;

        robot_store_->gt_v_s_(0) = msg->twist.twist.linear.x;
        robot_store_->gt_v_s_(1) = msg->twist.twist.linear.y;
        robot_store_->gt_v_s_(2) = msg->twist.twist.linear.z;

        // Test with GT orientation from Mocap
        // robot_store_->quaternion_.x() = msg->pose.pose.orientation.x;
        // robot_store_->quaternion_.y() = msg->pose.pose.orientation.y;
        // robot_store_->quaternion_.z() = msg->pose.pose.orientation.z;
        // robot_store_->quaternion_.w() = msg->pose.pose.orientation.w;

        // Quaterniond mocap_quaternion;
        gt_quaternion_.x() = msg->pose.pose.orientation.x;
        gt_quaternion_.y() = msg->pose.pose.orientation.y;
        gt_quaternion_.z() = msg->pose.pose.orientation.z;
        gt_quaternion_.w() = msg->pose.pose.orientation.w;
        EigenUtils::QuaternionToEuler(gt_quaternion_, gt_euler_);

        gt_p_ = robot_store_->gt_p_ - gt_p_offset_;
        gt_v_b_ = gt_quaternion_.normalized().toRotationMatrix().transpose() * robot_store_->gt_v_s_;
        gt_v_s_ = robot_store_->gt_v_s_;
    }

}
int main(int argc, char **argv)
{

    rclcpp::init(argc, argv);
    auto node = std::make_shared<robotSub::go1NewSub>("est_hardware_sub");
    rclcpp::spin(node);
    rclcpp::shutdown();

    return 0;
}
