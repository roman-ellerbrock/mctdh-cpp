find_package(UnitTest++ REQUIRED)
include_directories(${UTPP_INCLUDE_DIRS})

set(mctdh_tests
        test_hrep.cpp
        )

add_executable(test_mctdh ${mctdh_tests})
target_link_libraries(test_mctdh mctdh_lib)
target_link_libraries(test_mctdh UnitTest++)
target_compile_features(test_mctdh PRIVATE cxx_std_14)

# run unit tests as post build step
add_custom_command(TARGET test_mctdh
        POST_BUILD COMMAND test_mctdh
        COMMENT "Running unit tests")
