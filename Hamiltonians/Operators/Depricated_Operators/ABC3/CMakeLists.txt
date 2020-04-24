
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(ABC3_rst ABC3_rst.cpp ABC3_rst.f)
target_include_directories (ABC3_rst PUBLIC .)
target_link_libraries (ABC3_rst MCTDHlib)

