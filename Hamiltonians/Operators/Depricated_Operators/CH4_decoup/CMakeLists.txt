
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH4_decoup CH4_decoup.cpp CH4_decoup.f)
target_include_directories (CH4_decoup PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (CH4_decoup MCTDHlib)

