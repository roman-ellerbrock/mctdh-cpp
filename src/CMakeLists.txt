
add_library(mctdh_lib STATIC)
target_compile_features(mctdh_lib PRIVATE cxx_std_14)

target_compile_options(mctdh_lib PUBLIC "$<$<CONFIG:DEBUG>:${mctdh_DEBUG_FLAGS}>")
target_compile_options(mctdh_lib PRIVATE "$<$<CONFIG:RELEASE>:${mctdh_RELEASE_FLAGS}>")

target_sources(mctdh_lib PUBLIC
        Core/tCMFIntegrator.h
        Core/tCMFIntegrator.cpp
        Core/CMFIntegrator.h
        Core/CMFIntegrator.cpp
        Core/tEigenstates.h
        Core/tEigenstates.cpp
        Core/Eigenstates.h
        Core/Eigenstates.cpp

        Core/FortranOperator.h
        Core/FortranOperator.cpp
        Core/FortranSOP.h
        Core/FortranSOP.cpp

        Core/LayerInterface.h
        Core/LayerInterface.cpp
        Core/tLayerInterface.h
        Core/tLayerInterface.cpp
        Core/Hamiltonian.h
        Core/HamiltonianRepresentation.cpp
        Core/HamiltonianRepresentation.h
        Core/tHamiltonianRepresentation.cpp
        Core/tHamiltonianRepresentation.h
        Core/Wavefunction.h

        DVR/MatrixTensorTree.cpp
        DVR/MatrixTensorTree.h
        DVR/MatrixTensorTreeFunctions.h
        DVR/MatrixTensorTreeFunctions.cpp
        DVR/SymXMatrixTrees.cpp
        DVR/SymXMatrixTrees.h
        DVR/TreeGrids.h
        DVR/TDDVR.cpp
        DVR/TDDVR.h
        DVR/CDVR.cpp
        DVR/CDVR.h
        DVR/DeltaVTree.h
        DVR/DeltaVTree.cpp
        DVR/cdvr_functions.h
        DVR/cdvr_functions.cpp
#        DVR/GridRepresentation.cpp
#        DVR/GridRepresentation.h
        DVR/XMatrixTrees.cpp
        DVR/XMatrixTrees.h

        Parser/tyaml_parser.cpp
        Parser/tyaml_parser.h
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
