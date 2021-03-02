find_package(UnitTest++ REQUIRED)
include_directories(${UTPP_INCLUDE_DIRS})

set(mctdh_tests
        tests.cpp
        test_SL_Propagator.cpp
        test_ML_Propagator.cpp
        test_parser.cpp
        test_representation.cpp)

add_executable(tests ${mctdh_tests})

set(tests_DEBUG_FLAGS "-g -O3")
target_compile_options(tests PUBLIC "$<$<CONFIG:DEBUG>:${mctdh_tests_DEBUG_FLAGS}>")

target_link_libraries(tests mctdh_lib UnitTest++ QuTree::QuTree Hamiltonians yaml-cpp)

message(${PROJECT_SOURCE_DIR})

# run unit tests as post build step
add_custom_command(TARGET tests
        POST_BUILD COMMAND tests
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples
        COMMENT "Running unit tests")
