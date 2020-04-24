
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH4_radau CH4_radau_J0.cpp operatoren.ch4.radau.f)
target_include_directories (CH4_radau PUBLIC .)
target_link_libraries (CH4_radau MCTDHlib)

