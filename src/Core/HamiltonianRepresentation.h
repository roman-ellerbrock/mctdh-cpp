//
// Created by Roman Ellerbrock on 3/1/20.
//

#ifndef HAMILTONIANREPRESENTATION_H
#define HAMILTONIANREPRESENTATION_H
#include "Core/Hamiltonian.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"

class HamiltonianRepresentation {
public:
	HamiltonianRepresentation = default;
	~HamiltonianRepresentation = default;

	void build(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree) {
		Represent(hMats_, H, Psi, Psi, tree);
		Contraction(hContractions_, hMats_, Psi, Psi, tree);
	}

	SparseMatrixTree rho_;
	SparseMatrixTrees hMats_;
	SparseMatrixTrees hContractions_;

};


#endif //HAMILTONIANREPRESENTATION_H
