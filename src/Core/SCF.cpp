//
// Created by Hannes Hoppe on 16.06.22.
//

#include "SCF.h"
#include "TreeClasses/SpectralDecompositionTree.h"

namespace SCF {

    void SCF_Eigenstates(const Hamiltonian& H,
                         Tree& tree,
                         Wavefunction& Psi,
                         double eps,
                         size_t maxit){

        double error = 1e100;
        for(size_t i = 0; i < maxit and error > eps; ++i){
            SCF_Iteration(H, tree, Psi);
        }



    }

    void SCF_Iteration(const Hamiltonian& H,
                       Tree& tree,
                       Wavefunction& Psi){

        LinearizedSCFNodes ordered_nodes;
        ordered_nodes.linearize(tree);

        // there are always 2 nodes involved, except when applying the SCF to the top layer
        // normally, the current_node is the node the SCF should act on next, the last_node is
        // the previous node the SCF acted on. this last_node is the current top node, so
        // the first step for the SCF is to make the current node the new top node
        Node* current_node;
        Node* last_node{nullptr};

        // working wavefunction object


        // iterate over all nodes
        for(size_t i = 0; i < ordered_nodes.size(); ++i){
            current_node = ordered_nodes[i];
            // shift density matrix
            if(last_node != nullptr){
                SCF_ShiftTopNode(tree, Psi, last_node->address(), current_node->address());
            }

            // apply lanczos scheme

            std::swap(current_node, last_node);
        }
    }


    void SCF_ShiftTopNode(Tree& tree, Wavefunction& Psi, size_t last_node, size_t current_node){

        // QR decompose Tensor along given Dimension
        const int child_index = tree.nodes()[last_node].get().parentIdx();

        std::cout << child_index << std::endl;

        auto result = qr(Psi.attributes_[last_node], child_index);

        result.print(cout);
        exit(EXIT_SUCCESS);
    }

}

