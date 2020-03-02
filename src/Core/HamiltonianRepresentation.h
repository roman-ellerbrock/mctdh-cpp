//
// Created by Roman Ellerbrock on 3/1/20.
//

#ifndef HAMILTONIANREPRESENTATION_H
#define HAMILTONIANREPRESENTATION_H
#include "Core/Hamiltonian.h"
#include "Core/Wavefunction.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/MatrixTreeFunctions.h"

typedef SparseMatrixTrees<complex<double>> sMatrixTreeVector;

class HamiltonianRepresentation {
public:
	HamiltonianRepresentation(const Hamiltonian& H, Tree& tree): rho_(tree){
		for (const auto& M : H) {
			hMats_.emplace_back(SparseMatrixTreecd(M, tree));
			hContractions_.emplace_back(SparseMatrixTreecd(M, tree));
		}
	}

	~HamiltonianRepresentation() = default;

	void build(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree) {
		using namespace SparseMatrixTreeFunctions;
		using namespace MatrixTreeFunctions;

		Contraction(rho_, Psi, tree, true);
		Represent(hMats_, H, Psi, Psi, tree);
		Contraction(hContractions_, hMats_, Psi, Psi, tree);
	}

	void print(const Tree& tree, ostream& os = cout){
		os << "Rho:"<<endl;
		rho_.print(tree, os);
		os << "Matrix representations:" << endl;
		for (const auto& mat : hMats_) {
			mat.print(os);
		}
		for (const auto& con : hContractions_) {
			con.print(os);
		}
	}

	MatrixTreecd rho_;
	sMatrixTreeVector hMats_;
	sMatrixTreeVector hContractions_;

};


#endif //HAMILTONIANREPRESENTATION_H
