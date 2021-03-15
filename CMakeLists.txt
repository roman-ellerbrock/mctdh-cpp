###########################################################################
#
# --- mctdh++ ---
#
###########################################################################
cmake_minimum_required(VERSION 2.8)
project(mctdh-cpp C CXX)
enable_language(Fortran)
enable_language(C)
enable_language(CXX)
cmake_policy(SET "CMP0076" NEW)
cmake_policy(SET "CMP0040" NEW)
cmake_policy(SET "CMP0042" NEW)

# For debugging CMake
#set(CMAKE_VERBOSE_MAKEFILE ON)

include(vectorization.cmake)
set(mctdh_DEBUG_FLAGS "-g")
#set(FLAGS_VECTORIZE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect'")
#set(mctdh_RELEASE_FLAGS "${FLAGS_VECTORIZE} -O3 -ffast-math")
set(mctdh_RELEASE_FLAGS "-O3 -ffast-math")

if (APPLE)
    message("On OSX")
    set(mctdh_DEBUG_FLAGS "-mmacosx-version-min=10.14 ${mctdh_DEBUG_FLAGS}")
    set(mctdh_RELEASE_FLAGS "-mmacosx-version-min=10.14 ${mctdh_RELEASE_FLAGS}")
    set(CMAKE_MACOSX_RPATH 1)
elseif (UNIX)
    message("On UNIX")
else ()
    message("Other OS")
endif ()

separate_arguments(mctdh_DEBUG_FLAGS UNIX_COMMAND "${mctdh_DEBUG_FLAGS}")
separate_arguments(mctdh_RELEASE_FLAGS UNIX_COMMAND "${mctdh_RELEASE_FLAGS}")

find_package(QuTree REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories(${YAML_CPP_INCLUDE_DIR})
add_library(Yaml INTERFACE IMPORTED)

# build Hamiltonians
add_subdirectory("Hamiltonians")

# build mctdh_lib
add_subdirectory("src")

# build mctdh executable
add_executable(mctdh src/main.cpp)
target_compile_features(mctdh PRIVATE cxx_std_14)
target_compile_options(mctdh PUBLIC "$<$<CONFIG:DEBUG>:${mctdh_DEBUG_FLAGS}>")
target_compile_options(mctdh PRIVATE "$<$<CONFIG:RELEASE>:${mctdh_RELEASE_FLAGS}>")
target_include_directories(mctdh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(mctdh PUBLIC "-ffast-math" "-O3")
target_link_libraries(mctdh mctdh_lib QuTree::QuTree)

# build unit tests
add_subdirectory("tests")
