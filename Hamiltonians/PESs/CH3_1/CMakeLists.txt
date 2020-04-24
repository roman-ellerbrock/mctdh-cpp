
add_library(PESCH3_1 CH3Potential.cpp CH3Potential.f)
target_include_directories(PESCH3_1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PESCH3_1 MCTDHlib)

