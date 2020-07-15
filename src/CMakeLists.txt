
add_library(mctdh_lib SHARED)
target_compile_features(mctdh_lib PRIVATE cxx_std_14)

target_sources(mctdh_lib PUBLIC
        Core/CMFIntegrator.h
        Core/CMFIntegrator.cpp
        Core/Eigenstates.h
        Core/Eigenstates.cpp

        Core/FortranOperator.h
        Core/FortranOperator.cpp
        Core/FortranSOP.h
        Core/FortranSOP.cpp

        Core/LayerInterface.h
        Core/LayerInterface.cpp
        Core/Hamiltonian.h
        Core/HamiltonianRepresentation.cpp
        Core/HamiltonianRepresentation.h
        Core/Wavefunction.h

        DVR/ExplicitEdgeWavefunction.cpp
        DVR/ExplicitEdgeWavefunction.h
        DVR/TreeGrids.h
        DVR/TDDVR.cpp
        DVR/TDDVR.h
        DVR/CDVR.cpp
        DVR/CDVR.h
        DVR/DeltaVTree.h
        DVR/DeltaVTree.cpp
        DVR/cdvr_functions.h
        DVR/cdvr_functions.cpp
        DVR/XMatrixTrees.cpp
        DVR/XMatrixTrees.h

        Parser/yaml_parser.cpp
        Parser/yaml_parser.h
        )

target_include_directories(mctdh_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )


target_link_libraries(mctdh_lib
        Hamiltonians
        yaml-cpp
        QuTree::QuTree)
