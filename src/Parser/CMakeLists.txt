


add_library(parser)

target_sources(parser
        PUBLIC
        yaml_parser.cpp)

target_include_directories(parser
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )


target_compile_features(parser PRIVATE cxx_std_14)
