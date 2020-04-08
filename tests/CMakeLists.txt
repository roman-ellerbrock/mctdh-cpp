find_package(UnitTest++ REQUIRED)
include_directories(${UTPP_INCLUDE_DIRS})

set(mctdh_tests
        test_CDVR.cpp
        test_ExplicitEdgeWavefunction.cpp
        test_TDDVR.cpp
        test_hrep.cpp
        test_Integrator.cpp
        tests.cpp
        )

add_executable(tests ${mctdh_tests})
target_link_libraries(tests mctdh_lib UnitTest++ QuTree::QuTree)

# run unit tests as post build step
add_custom_command(TARGET tests
        POST_BUILD COMMAND tests
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples
        COMMENT "Running unit tests")
