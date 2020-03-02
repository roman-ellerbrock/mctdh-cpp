//
// Created by Roman Ellerbrock on 3/1/20.
//

#ifndef HAMILTONIANREPRESENTATION_H
#define HAMILTONIANREPRESENTATION_H
#include "Core/Hamiltonian.h"
#include "Core/Wavefunction.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SpectralDecompositionTree.h"

typedef SparseMatrixTrees<complex<double>> sMatrixTreeVector;

class HamiltonianRepresentation {
public:
	HamiltonianRepresentation(const Hamiltonian& H, Tree& tree): rho_(tree),
	rho_decomposition_(tree), rho_inverse_(tree){
		for (const auto& M : H) {
			hMats_.emplace_back(SparseMatrixTreecd(M, tree));
			hContractions_.emplace_back(SparseMatrixTreecd(M, tree));
		}
	}

	~HamiltonianRepresentation() = default;

	void build(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree) {
		using namespace SparseMatrixTreeFunctions;
		using namespace MatrixTreeFunctions;

		/// Calculate density matrix tree
		Contraction(rho_, Psi, tree, true);

		/// Density matrix tree decomposition
		rho_decomposition_.Calculate(rho_, tree);

		/// Calculate inverse density matrix tree
		rho_inverse_ = rho_decomposition_.Invert(tree);

		/// Calculate h-matrix trees
		Represent(hMats_, H, Psi, Psi, tree);

		/// Calculate h-matrix tree contractions
		Contraction(hContractions_, hMats_, Psi, Psi, tree);
	}

	void print(const Tree& tree, ostream& os = cout){
		os << "Rho:"<<endl;
		rho_.print(tree, os);
		os << "Matrix representations:" << endl;
		for (const auto& mat : hMats_) {
			mat.print(os);
		}
		os << "Matrix Contractions:" << endl;
		for (const auto& con : hContractions_) {
			con.print(os);
		}
	}

	MatrixTreecd rho_;
	SpectralDecompositionTreecd rho_decomposition_;
	MatrixTreecd rho_inverse_;
	sMatrixTreeVector hMats_;
	sMatrixTreeVector hContractions_;

};

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node) {
	using namespace SparseMatrixTreeFunctions;

	dPhi.Zero();
	// For each part in the hamiltonian
	for (size_t l = 0; l < H.size(); l++) {
		// Gather information on what is active at this node.
		// That is the MPO and the CDVR-Matrices acting on this node
		const MLOcd& M = H[l];
		// d is the number of CDVR children at this node
		// Apply the hmatrices for this part
		const auto& hmat = hRep.hMats_[l];
		// Proceed only if H-matrices or CDVR acts on this node
		if (!hmat.Active(node)) { continue; }

		// Get the factor of the part
		const complex<double> coeff = H.Coeff(l);
		// Copy the Acoeffs and multyply with the coefficient
		Tensorcd Psi(Phi, coeff);
		Psi = SparseMatrixTreeFunctions::Apply(hmat, Psi, M, node);

		// Hole-Matrix
		// Multiply with H-Hole-matrix at this layer_
		const auto& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.Active(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
			dPhi += Psi;
		}
	}

	// Now the part in the EOM that only has to be done once:
	// Inverse Densitymatrix and (1-P) projector for SPF-type EOM
	if (!node.isToplayer()) {
		dPhi = ProjectOut(dPhi, Phi);

		// Multiply with inverse single-particle density matrix
		const auto& rhoinv = hRep.rho_inverse_[node];
		dPhi = multStateAB(rhoinv, dPhi);
	}
}

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree) {
	hRep.build(H, Psi, tree);
	for (const Node& node : tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
}

#endif //HAMILTONIANREPRESENTATION_H
