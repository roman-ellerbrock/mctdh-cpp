
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(GenericKEO GenericKEO.cpp)
target_include_directories (GenericKEO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (GenericKEO MCTDHlib)

