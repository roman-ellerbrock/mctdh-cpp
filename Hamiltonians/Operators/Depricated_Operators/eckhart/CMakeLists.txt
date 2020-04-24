
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(eckhart eckhart.cpp eckhart.f)
target_include_directories (eckhart PUBLIC .)
target_link_libraries (eckhart MCTDHlib)

