
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(Hubbard hamiltonianHubbardSQR.cpp)
target_include_directories (Hubbard PUBLIC .)
target_link_libraries (Hubbard MCTDHlib)

