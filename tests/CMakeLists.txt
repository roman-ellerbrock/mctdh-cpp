find_package(UnitTest++ REQUIRED)
include_directories(${UTPP_INCLUDE_DIRS})

set(mctdh_tests
        test_hrep.cpp
        test_Integrator.cpp
        tests.cpp
        )

add_executable(tests ${mctdh_tests})
target_link_libraries(tests mctdh_lib UnitTest++)
target_compile_features(tests PRIVATE cxx_std_14)

# run unit tests as post build step
add_custom_command(TARGET test
        POST_BUILD COMMAND tests
        COMMENT "Running unit tests")
