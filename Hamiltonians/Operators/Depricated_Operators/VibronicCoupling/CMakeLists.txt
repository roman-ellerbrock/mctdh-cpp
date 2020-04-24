
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(VibronicCoupling hamiltonianCnv.cpp)
target_include_directories (VibronicCoupling PUBLIC .)
target_link_libraries (VibronicCoupling MCTDHlib)

