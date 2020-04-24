
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(HCH4_rst HCH4_rst.cpp HCH4_rst.f)
target_include_directories (HCH4_rst PUBLIC .)
target_link_libraries (HCH4_rst MCTDHlib)

