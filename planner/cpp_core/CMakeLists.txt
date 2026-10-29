cmake_minimum_required(VERSION 3.16)
project(pingpong_planner_core VERSION 0.1.0 LANGUAGES CXX)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(pingpong_planner_core
  src/trajectory_predictor.cpp
  src/racket_target_planner.cpp
)
target_compile_features(pingpong_planner_core PUBLIC cxx_std_17)
set_target_properties(pingpong_planner_core PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(pingpong_planner_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(pingpong_planner_core PUBLIC Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(pingpong_planner_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
  find_package(GTest REQUIRED)
  add_executable(test_trajectory_and_stage3 tests/test_trajectory_and_stage3.cpp)
  target_link_libraries(test_trajectory_and_stage3 PRIVATE
    pingpong_planner_core GTest::gtest_main
  )
  include(GoogleTest)
  gtest_discover_tests(test_trajectory_and_stage3)
endif()
