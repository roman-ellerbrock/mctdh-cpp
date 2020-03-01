

target_sources(mctdh
        PUBLIC
        main.cpp
        Core/Hamiltonian.h
        Core/Wavefunction.h
        Parser/yaml_parser.cpp
        Parser/yaml_parser.h
        )


add_subdirectory("Parser")


