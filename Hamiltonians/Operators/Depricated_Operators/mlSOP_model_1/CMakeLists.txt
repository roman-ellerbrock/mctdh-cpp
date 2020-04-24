
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(mlSOPm1 mlSOP_model_1.cpp)
target_include_directories (mlSOPm1 PUBLIC .)
target_link_libraries (mlSOPm1 MCTDHlib)

