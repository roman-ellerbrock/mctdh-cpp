//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "Core/tHamiltonianRepresentation.h"

template<typename T>
Tensor<T> Apply(const SOP<T>& H, const Tensor<T>& Phi,
	const tHamiltonianRepresentation<T>& hRep,
	const Node& node) {

	Tensor<T> dPhi(Phi.shape());
	Tensor<T> Psi(Phi.shape());
	for (size_t l = 0; l < H.size(); l++) {
		const auto& hmat = hRep.hMats_[l];
		if (!hmat.isActive(node)) { continue; }
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			Psi[i] = H.coeff(l) * Phi[i];
		}

//		Psi = H.Coeff(l) * Phi;
//		Tensorcd Psi(Phi, H.Coeff(l));
		Psi = TreeFunctions::apply(hmat, Psi, H[l], node);

		// Multiply with hole-matrix
		const auto& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.isActive(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
			dPhi += Psi;
		}
	}

	if (H.hasV) {
		dPhi += hRep.cdvr_.Apply(Phi, hRep.rho_decomposition_[node], node);
	}

	return dPhi;
}

template<typename T>
Matrix<T> Expectation(const tHamiltonianRepresentation<T>& hRep,
	const TensorTree<T>& Psi, const SOP<T>& H, const Tree& tree) {

	const Node& top = tree.topNode();
	auto HPhi = Apply(H, Psi[top], hRep, top);
	return Psi[top].dotProduct(HPhi);
}

template<typename T>
void LayerDerivative(Tensor<T>& dPhi, double time, const Tensor<T>& Phi,
	const SOP<T>& H, const tHamiltonianRepresentation<T>& hRep,
	const Node& node, T propagation_phase) {

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

	dPhi *= propagation_phase;
}

template<typename T>
void Derivative(TensorTree<T>& dPsi, tHamiltonianRepresentation<T>& hRep,
	double time, const TensorTree<T>& Psi, const SOP<T>& H,
	const Tree& tree, T propagation_phase) {
	hRep.build(H, Psi, tree);
	for (const Node& node : tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
}

template class tHamiltonianRepresentation<double>;

