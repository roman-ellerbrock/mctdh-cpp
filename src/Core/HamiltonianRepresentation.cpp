//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "Core/HamiltonianRepresentation.h"

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep,
	const Node& node, double time) {

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
	double time) {

	const Node& top = tree.TopNode();
	auto dPhi = Apply(H, Psi[top], hRep, top, time);
	return Psi[top].DotProduct(dPhi);
}

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase) {

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
}

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase) {
	hRep.build(H, Psi, tree);
	for (const Node& node : tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
	dPsi *= -propagation_phase * QM::im;
}

