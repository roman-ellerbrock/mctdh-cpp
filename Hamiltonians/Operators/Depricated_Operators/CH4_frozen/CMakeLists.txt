
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH4_frozen CH4_frozen.cpp CH4_frozen.f)
target_include_directories (CH4_frozen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (CH4_frozen MCTDHlib)

