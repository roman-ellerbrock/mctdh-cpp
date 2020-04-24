
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CDVRModel3 CDVRModel3.cpp)
target_include_directories (CDVRModel3 PUBLIC .)
target_link_libraries (CDVRModel3 MCTDHlib)

