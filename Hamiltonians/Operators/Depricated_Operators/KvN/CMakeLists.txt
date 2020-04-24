
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(KvN hamiltonKvN.cpp)
target_include_directories (KvN PUBLIC .)
target_link_libraries (KvN MCTDHlib)

