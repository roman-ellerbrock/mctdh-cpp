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

#set(CMAKE_VERBOSE_MAKEFILE ON)

#find_package(qutree)

find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})
add_library(Yaml INTERFACE IMPORTED)

find_package(Eigen3 REQUIRED)
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("QDlib/include")

add_subdirectory("Hamiltonians")

# Library
add_executable(mctdh src/main.cpp)
target_compile_features(mctdh PRIVATE cxx_std_14)
target_include_directories(mctdh
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/QDlib/include
        )
target_link_directories(mctdh
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/QDlib/cmake-build-debug
        )
add_subdirectory("src")
target_link_libraries(mctdh
        QDlib
        Hamiltonians
        yaml-cpp)

# Executable
#add_executable(mctdh src/main.cpp)
#target_include_directories(mctdh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
#target_link_libraries(mctdh mctdh_lib)

#add_subdirectory("tests")
