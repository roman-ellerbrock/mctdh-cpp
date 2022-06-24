//
// Created by Hannes Hoppe on 16.06.22.
//

#ifndef MCTDH_CPP_SCF_H
#define MCTDH_CPP_SCF_H

#include "Core/Hamiltonian.h"
#include "Core/Wavefunction.h"
#include "TreeShape/LinearizedSCFNodes.h"

// this defines an SCF procedure to find eigenstates of a given hamiltonian

namespace SCF {

    void SCF_Eigenstates(const Hamiltonian &H,
                         const Tree &tree,
                         Wavefunction &Psi,
                         double eps = 1e-12,
                         size_t maxit = 100);

    void SCF_Iteration(const Hamiltonian &H,
                       Tree &tree,
                       Wavefunction &Psi_In);

    void SCF_ShiftTopNode(Tree& tree,
                          Wavefunction& Psi,
                          size_t from_node,
                          size_t to_node);

}
#endif //MCTDH_CPP_SCF_H
