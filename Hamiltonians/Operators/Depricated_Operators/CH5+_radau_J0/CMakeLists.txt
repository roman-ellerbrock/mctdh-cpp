
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH5+_radau CH5+_radau_J0.cpp operatoren.ch5+.f)
target_include_directories (CH5+_radau PUBLIC .)
target_link_libraries (CH5+_radau MCTDHlib)

