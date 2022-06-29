//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "Core/HamiltonianRepresentation.h"

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep,
	const Node& node) {

	Tensorcd dPhi(Phi.shape());
	Tensorcd Psi(Phi.shape());
	for (size_t l = 0; l < H.size(); l++) {
		const auto& hmat = hRep.hMats_[l];
		if (!hmat.isActive(node)) { continue; }
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			Psi[i] = H.coeff(l) * Phi[i];
		}

		Psi = TreeFunctions::apply(hmat, Psi, H[l], node);

		// Multiply with hole-matrix
		const SparseMatrixTreecd& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.isActive(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
			dPhi += Psi;
		}
	}

	if (H.hasV) {
		dPhi += hRep.cdvr_.apply(Phi, hRep.rho_decomposition_[node], node);
	}

	return dPhi;
}

Matrixcd Expectation(const HamiltonianRepresentation& hRep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree) {

	const Node& top = tree.topNode();
	auto HPhi = Apply(H, Psi[top], hRep, top);
	return Psi[top].dotProduct(HPhi);
}

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase) {

	// dPhi = -i*rho^-1 * (1-P) * <Phi|H|Phi>
	// Apply Hamiltonian
	dPhi = Apply(H, Phi, hRep, node);

	// Inverse Densitymatrix and (1-P) projector for SPF-type EOM
	if (!node.isToplayer()) {
		// (1-P) projector
		dPhi = projectOut(dPhi, Phi);

		// Multiply with inverse single-particle density matrix
		const auto& rhoinv = hRep.rho_inverse_[node];
		dPhi = multStateAB(rhoinv, dPhi);
	}

	dPhi *= -propagation_phase * QM::im;
}

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase) {
	hRep.build(H, Psi, tree, time);
	for (const Node& node : tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
}

////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////

void symLayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase) {

//	dPhi = TreeFunctions::symApply(Phi, hRep.hMatSets_, H, node);
	dPhi *= -propagation_phase * QM::im;
}

void symDerivative(MatrixTensorTree& dPsi, HamiltonianRepresentation& hRep,
	double time, const MatrixTensorTree& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase) {

	for (const Node& node : tree) {
		symLayerDerivative(dPsi.first[node], time, Psi.first[node],
			H, hRep, node, propagation_phase);
	}

}
