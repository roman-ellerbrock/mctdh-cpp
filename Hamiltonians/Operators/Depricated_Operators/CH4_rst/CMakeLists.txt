
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH4_rst CH4_rst.cpp CH4_rst.f)
target_include_directories (CH4_rst PUBLIC .)
target_link_libraries (CH4_rst MCTDHlib)

