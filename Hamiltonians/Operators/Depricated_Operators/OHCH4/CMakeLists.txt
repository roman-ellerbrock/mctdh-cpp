
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(OHCH4 OHCH4.cpp OHCH4.f)
target_include_directories (OHCH4 PUBLIC .)
target_link_libraries (OHCH4 MCTDHlib)

