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
#include "Util/QMConstants.h"

typedef SparseMatrixTrees<complex<double>> sMatrixTreeVector;

class HamiltonianRepresentation {
public:
	HamiltonianRepresentation(const Hamiltonian& H, const Tree& tree)
		: rho_(tree),
		  rho_decomposition_(tree), rho_inverse_(tree) {
		for (const auto& M : H) {
			hMats_.emplace_back(SparseMatrixTreecd(M, tree));
			hContractions_.emplace_back(SparseMatrixTreecd(M, tree));
		}
	}

	~HamiltonianRepresentation() = default;

	void build(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree) {
		/// Calculate density matrix tree
		MatrixTreeFunctions::Contraction(rho_, Psi, tree, true);

		/// Density matrix tree decomposition
		rho_decomposition_.Calculate(rho_, tree);

		/// Calculate inverse density matrix tree
		rho_inverse_ = rho_decomposition_.Invert(tree);

		/// Calculate h-matrix trees
		SparseMatrixTreeFunctions::Represent(hMats_, H, Psi, Psi, tree);

		/// Calculate h-matrix tree contractions
		SparseMatrixTreeFunctions::Contraction(hContractions_, hMats_, Psi, Psi, tree);
	}

	void print(const Tree& tree, ostream& os = cout) {
		os << "Rho:" << endl;
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

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep,
	const Node& node, double time = 0.) {

	Tensorcd dPhi(Phi.shape());
	for (size_t l = 0; l < H.size(); l++) {
		const auto& hmat = hRep.hMats_[l];
		if (!hmat.Active(node)) { continue; }

		Tensorcd Psi(Phi, H.Coeff(l));
		Psi = SparseMatrixTreeFunctions::Apply(hmat, Psi, H[l], node);

		// Multiply with hole-matrix
		const auto& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.Active(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
			dPhi += Psi;
		}
	}
	return dPhi;
}

Matrixcd Expectation(const HamiltonianRepresentation& hRep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree,
	double time = 0.) {

	const Node& top = tree.TopNode();
	auto dPhi = Apply(H, Psi[top], hRep, top, time);
	return Psi[top].DotProduct(dPhi);
}

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase = 1.) {

	// Apply Hamiltonian
	dPhi = Apply(H, Phi, hRep, node, time);

	// Inverse Densitymatrix and (1-P) projector for SPF-type EOM
	if (!node.isToplayer()) {
		// (1-P) projector
		dPhi = ProjectOut(dPhi, Phi);

		// Multiply with inverse single-particle density matrix
		const auto& rhoinv = hRep.rho_inverse_[node];
		dPhi = multStateAB(rhoinv, dPhi);
	}

	dPhi *= propagation_phase * QM::im;
}

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase = 1.) {
	hRep.build(H, Psi, tree);
	for (const Node& node : tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
}


#endif //HAMILTONIANREPRESENTATION_H
