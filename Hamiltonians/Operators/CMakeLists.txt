
#add_library(Operators CWrapper.f Operators.cpp GenericKEO/GenericKEO.cpp GenericKEO/GenericKEO.h Operator_lib.cpp Operator_lib.h NumberPartitioning/NumberPartitioning.cpp NumberPartitioning/NumberPartitioning.h IntegerFactorisation/IntegerFactorisation.h IntegerFactorisation/IntegerFactorisation.cpp UnsupervisedLearning/PictureInterface.cpp UnsupervisedLearning/PictureInterface.h)
#target_include_directories(Operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#############################################
# List all available Operators here
# Note: Remember to link them to the Operators-library
#############################################
add_library(Operators Operators.h
        CoupledHO/CoupledHO.cpp
        CoupledHO/CoupledHO.h
        )

target_compile_features(Operators PRIVATE cxx_std_14)
target_link_libraries(Operators QuTree::QuTree)

#target_link_libraries(Operators LINK_PUBLIC
#        CoupledHO
#        )

#add_subdirectory("CH3_meanfield")
#add_subdirectory("CH3_quasiexact")
#add_subdirectory("CoupledHO")
#add_subdirectory("GateOperators")
#add_subdirectory("Ising")
#add_subdirectory("StandardKEO")
#add_subdirectory("NumberPartitioning")
#add_subdirectory("IntegerFactorisation")
#add_subdirectory("scattering")
#add_subdirectory("UnsupervisedLearning")
#
#target_link_libraries(Operators LINK_PUBLIC
#        CH3_quasiexact
#        CH3_meanfield
#        CoupledHO
#        GateOperators
#        Ising
#        NOCl
#        StandardKEO
#        NumberPartitioning
#        IntegerFactorisation
#        scattering
#        UnsupervisedLearning
#        )