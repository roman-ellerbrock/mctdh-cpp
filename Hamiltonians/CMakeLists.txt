
# Compile the Operators
add_subdirectory("Operators")

# Compile the Potential Energy Surfaces
add_subdirectory("PESs")


# Link the libraries "Operators" and "PESs" into a single library "Hamiltonians"
add_library(Hamiltonians Hamiltonians.cpp Hamiltonians.h)
target_compile_features(Hamiltonians PRIVATE cxx_std_14)
target_link_libraries(Hamiltonians LINK_PUBLIC Operators PESs)
target_link_libraries(Hamiltonians LINK_PUBLIC QuTree::QuTree)
target_include_directories(Hamiltonians PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

