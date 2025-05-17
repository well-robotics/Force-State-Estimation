cmake_minimum_required(VERSION 3.10)

project(go1_example)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(std_msgs REQUIRED)
find_package(nav_msgs REQUIRED)


find_package(decentral_legged_est REQUIRED)

# simulation dependencies
find_package(unitree_go REQUIRED)
find_package(unitree_api REQUIRED)

# hardware dependencies
find_package(sensor_msgs REQUIRED)
find_package(optitrack_broadcast REQUIRED)

include_directories(SYSTEM
  ${EIGEN3_INCLUDE_DIR}
  ${decentral_legged_est_INCLUDE_DIRS}
  include
)
# aux_source_directory(. DIR_SRCS)

add_library(code_gen
  src/Expressions/u_map_b1_description.hh
  src/Expressions/u_map_b1_description.cc
  src/Expressions/Ge_vec_b1_description.hh
  src/Expressions/Ge_vec_b1_description.cc
  src/Expressions/Mmat1_b1_description.hh
  src/Expressions/Mmat1_b1_description.cc
  src/Expressions/Mmat2_b1_description.hh
  src/Expressions/Mmat2_b1_description.cc
  src/Expressions/Mmat3_b1_description.hh
  src/Expressions/Mmat3_b1_description.cc
  src/Expressions/Mmat4_b1_description.hh
  src/Expressions/Mmat4_b1_description.cc
  src/Expressions/Mmat5_b1_description.hh
  src/Expressions/Mmat5_b1_description.cc
  src/Expressions/Mmat6_b1_description.hh
  src/Expressions/Mmat6_b1_description.cc
  src/Expressions/Mmat7_b1_description.hh
  src/Expressions/Mmat7_b1_description.cc
  src/Expressions/Mmat8_b1_description.hh
  src/Expressions/Mmat8_b1_description.cc
  src/Expressions/Mmat9_b1_description.hh
  src/Expressions/Mmat9_b1_description.cc
  src/Expressions/Mmat10_b1_description.hh
  src/Expressions/Mmat10_b1_description.cc
  src/Expressions/Mmat11_b1_description.hh
  src/Expressions/Mmat11_b1_description.cc
  src/Expressions/Mmat12_b1_description.hh
  src/Expressions/Mmat12_b1_description.cc
  src/Expressions/Mmat13_b1_description.hh
  src/Expressions/Mmat13_b1_description.cc
  src/Expressions/Mmat14_b1_description.hh
  src/Expressions/Mmat14_b1_description.cc
  src/Expressions/C_mat_go1.hh
  src/Expressions/C_mat_go1.cc
  src/Expressions/FL_foot.cc
  src/Expressions/FR_foot.cc
  src/Expressions/RL_foot.cc
  src/Expressions/RR_foot.cc
  src/Expressions/J_FL.cc
  src/Expressions/J_FR.cc
  src/Expressions/J_RL.cc
  src/Expressions/J_RR.cc
  src/Expressions/FL_foot.hh
  src/Expressions/FR_foot.hh
  src/Expressions/RL_foot.hh
  src/Expressions/RR_foot.hh
  src/Expressions/J_FL.hh
  src/Expressions/J_FR.hh
  src/Expressions/J_RL.hh
  src/Expressions/J_RR.hh
)

# realtime_est

add_executable(grf_sub src/go1NewSub.cpp)
target_link_libraries(grf_sub code_gen OsqpEigen::OsqpEigen)
ament_target_dependencies(grf_sub rclcpp decentral_legged_est std_msgs sensor_msgs nav_msgs optitrack_broadcast)

install(TARGETS grf_sub DESTINATION lib/${PROJECT_NAME})

install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}/
)
install(FILES config/parameters_go1_new.yaml
  DESTINATION share/${PROJECT_NAME}/config)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()

