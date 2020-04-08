###########################################################################
#
# --- mctdh++ ---
#
###########################################################################
cmake_minimum_required(VERSION 2.8)
project(mctdh-cpp C CXX)
enable_language(C)
enable_language(CXX)
cmake_policy(SET "CMP0076" NEW)
cmake_policy(SET "CMP0040" NEW)

# For debugging CMake
#set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(QuTree REQUIRED)

find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})
add_library(Yaml INTERFACE IMPORTED)

#find_package(Eigen3 REQUIRED)
#include_directories("${EIGEN3_INCLUDE_DIR}")

#include_directories("QDlib/include")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/QDlib/cmake-build-debug")

# build Hamiltonians
add_subdirectory("Hamiltonians")

# build mctdh_lib
add_subdirectory("src")

# build mctdh executable
add_executable(mctdh src/main.cpp)
target_compile_features(mctdh PRIVATE cxx_std_14)
target_include_directories(mctdh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(mctdh mctdh_lib QuTree::QuTree)

# build unit tests
add_subdirectory("tests")
