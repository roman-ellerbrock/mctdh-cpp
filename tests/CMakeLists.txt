include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/heads/main.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()


set(mctdh_tests
        tests.cpp
        test_TDDVR.cpp
        )
add_executable(tests ${mctdh_tests})

set(tests_DEBUG_FLAGS "-g -O3")
target_compile_options(tests PUBLIC "$<$<CONFIG:DEBUG>:${mctdh_tests_DEBUG_FLAGS}>")

target_link_libraries(tests mctdh_lib gtest_main QuTree::QuTree Hamiltonians yaml-cpp)

message(${PROJECT_SOURCE_DIR})

# run unit tests as post build step
add_custom_command(TARGET tests
        POST_BUILD COMMAND tests
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples
        COMMENT "Running unit tests")
