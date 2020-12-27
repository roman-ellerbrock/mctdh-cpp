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
		if (!hmat.Active(node)) { continue; }
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			Psi[i] = H.Coeff(l) * Phi[i];
		}

//		Psi = H.Coeff(l) * Phi;
//		Tensorcd Psi(Phi, H.Coeff(l));
		Psi = TreeFunctions::Apply(hmat, Psi, H[l], node);

		// Multiply with hole-matrix
		const auto& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.Active(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
//			for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
//				dPhi[i] += Psi[i];
//			}
			dPhi += Psi;
		}
	}

	if (H.hasV) {
		auto x_sqrt_rho = sqrt(hRep.rho_decomposition_[node]);
		x_sqrt_rho.second = Regularize(x_sqrt_rho.second, 1e-6);
		auto sqrt_rho = toMatrix(x_sqrt_rho);
		auto VPhi = hRep.cdvr_.Apply(Phi, sqrt_rho, node);
		dPhi += VPhi;
	}

	return dPhi;
}

Matrixcd Expectation(const HamiltonianRepresentation& hRep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree) {

	const Node& top = tree.TopNode();
	auto HPhi = Apply(H, Psi[top], hRep, top);
	return Psi[top].DotProduct(HPhi);
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
		dPhi = ProjectOut(dPhi, Phi);

		// Multiply with inverse single-particle density matrix
		const auto& rhoinv = hRep.rho_inverse_[node];
		dPhi = multStateAB(rhoinv, dPhi);
	}

	dPhi *= -propagation_phase * QM::im;
}

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase) {
	hRep.build(H, Psi, tree);
	for (const Node& node : tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
}

