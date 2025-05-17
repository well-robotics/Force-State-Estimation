cmake_minimum_required(VERSION 3.10)

project(decentral_legged_est)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(custom_msgs REQUIRED)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src/Expressions include)

add_library(codegen_dyn src/Expressions/u_map_b1_description.hh
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
                        src/Expressions/J_FR.hh
                        src/Expressions/J_FR.cc
                        src/Expressions/J_FL.hh
                        src/Expressions/J_FL.cc
                        src/Expressions/J_RR.hh
                        src/Expressions/J_RR.cc
                        src/Expressions/J_RL.hh
                        src/Expressions/J_RL.cc
)

add_library(${PROJECT_NAME} SHARED 
src/Spline/Bezier_simple.cpp
src/MheSrb.cpp
src/DecentralEst.cpp
src/EstSub.cpp)

target_link_libraries(${PROJECT_NAME} OsqpEigen::OsqpEigen codegen_dyn)
ament_target_dependencies(${PROJECT_NAME} rclcpp custom_msgs sensor_msgs geometry_msgs) 

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS ${PROJECT_NAME}
  EXPORT export_${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION include/${PROJECT_NAME}/
)

ament_export_include_directories(include)
ament_export_libraries(${PROJECT_NAME})
ament_export_dependencies(rclcpp custom_msgs sensor_msgs geometry_msgs)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
#MPCExample
# add_executable(mhe src/mheSRB.cpp)
# target_link_libraries(mhe OsqpEigen::OsqpEigen transforms forward_kinematics)
