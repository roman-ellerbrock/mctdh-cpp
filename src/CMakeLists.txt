
add_library(mctdh_lib SHARED)
target_compile_features(mctdh_lib PRIVATE cxx_std_14)

target_sources(mctdh_lib
        PUBLIC
        Core/CMFIntegrator.h
        Core/CMFIntegrator.cpp
        Core/LayerInterface.h
        Core/LayerInterface.cpp
        Core/Hamiltonian.h
        Core/HamiltonianRepresentation.cpp
        Core/HamiltonianRepresentation.h
        Core/Wavefunction.h
        Parser/yaml_parser.cpp
        Parser/yaml_parser.h
        )


target_include_directories(mctdh_lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )


target_link_libraries(mctdh_lib
        QDlib
        Hamiltonians
        yaml-cpp)
