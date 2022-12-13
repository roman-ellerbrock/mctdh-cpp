
add_library(mctdh_lib STATIC)
target_compile_features(mctdh_lib PRIVATE cxx_std_14)

target_compile_options(mctdh_lib PUBLIC "$<$<CONFIG:DEBUG>:${mctdh_DEBUG_FLAGS}>")
target_compile_options(mctdh_lib PRIVATE "$<$<CONFIG:RELEASE>:${mctdh_RELEASE_FLAGS}>")

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
        Core/QuadraticSOP.cpp
        Core/QuadraticSOP.h
        Core/QuadraticSVD.cpp
        Core/QuadraticSVD.h
        Core/Wavefunction.h
        Core/Wavefunction.cpp
        Core/SCF.cpp
        Core/SCF.h

        ../../QuTree/src/DVR/MatrixTensorTree.cpp
        ../../QuTree/include/DVR/MatrixTensorTree.h
        ../../QuTree/include/DVR/MatrixTensorTreeFunctions.h
        ../../QuTree/src/DVR/MatrixTensorTreeFunctions.cpp
        ../../QuTree/src/DVR/SymXMatrixTrees.cpp
        ../../QuTree/include/DVR/SymXMatrixTrees.h
        ../../QuTree/include/DVR/TreeGrids.h
        ../../QuTree/src/DVR/TDDVR.cpp
        ../../QuTree/include/DVR/TDDVR.h
        ../../QuTree/src/DVR/CDVR.cpp
        ../../QuTree/include/DVR/CDVR.h
        ../../QuTree/include/DVR/DeltaVTree.h
        ../../QuTree/src/DVR/DeltaVTree.cpp
        ../../QuTree/include/DVR/cdvr_functions.h
        ../../QuTree/src/DVR/cdvr_functions.cpp
#        DVR/GridRepresentation.cpp
#        DVR/GridRepresentation.h
        ../../QuTree/src/DVR/XMatrixTrees.cpp
        ../../QuTree/include/DVR/XMatrixTrees.h

        Parser/yaml_parser.cpp
        Parser/yaml_parser.h

        Utility/Overlaps.cpp
        Utility/normal_modes.cpp
        Utility/NewtonRaphson.cpp
        )

target_include_directories(mctdh_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

target_link_libraries(mctdh_lib
        Hamiltonians
        yaml-cpp
        QuTree::QuTree)
