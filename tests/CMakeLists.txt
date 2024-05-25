cmake_minimum_required(VERSION 2.6)
 
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN_INCLUDE_DIR})

find_package(Threads)

# Link runTests with what we want to test and the GTest and pthread library
add_executable(run_tests.x test_inv.cpp run_tests.cpp)
target_link_libraries(run_tests.x ${GTEST_LIBRARIES} pthread Eigen3::Eigen)