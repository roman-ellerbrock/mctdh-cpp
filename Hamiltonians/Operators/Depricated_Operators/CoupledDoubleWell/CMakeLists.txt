
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CoupledDoubleWell CoupledDoubleWell.cpp CoupledDoubleWell.h)
target_include_directories (CoupledDoubleWell PUBLIC .)
target_link_libraries (CoupledDoubleWell MCTDHlib)

