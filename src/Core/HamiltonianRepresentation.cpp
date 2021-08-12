//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "Core/HamiltonianRepresentation.h"

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep,
	const Node& node) {

	Tensorcd dPhi(Phi.shape());
	Tensorcd Psi(Phi.shape());
	/// For every term in the Hamiltonian
	for (size_t l = 0; l < H.size(); l++) {
		const auto& hmat = hRep.hMats_[l];
		if (!hmat.isActive(node)) { continue; }

		/// Initialize Psi via |Psi> = c_l * |Phi>
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			Psi[i] = H.coeff(l) * Phi[i];
		}

		/// Apply matrix representation of H for every term l
		Psi = TreeFunctions::apply(hmat, Psi, H[l], node);

		/// Multiply with Mean-field matrix
		const auto& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.isActive(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
			dPhi += Psi;
		}
	}

	/// Add CDVR if potential energy surface is used
	if (H.hasV) {
		dPhi += hRep.cdvr_.Apply(Phi, hRep.rho_decomposition_[node], node);
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
	/// \brief This routine evaluates derivative d/dt |Phi^p> = -i*rho^-1 * (1-P) * <Psi(p)|H|Psi(p)>|Phi^p>

	/// 1.) Apply Hamiltonian, a.) H-matrices and b.) H-mean-field matrices
	dPhi = Apply(H, Phi, hRep, node);

	if (!node.isToplayer()) {
		/// 2.) Project out rotation in active space, |hPhi^p > = (1-P)|hPhi^p>
		dPhi = projectOut(dPhi, Phi);

		/// 3.) Multiply by inverse density matrix rho^{(p),-1}|hPhi^p>
		const auto& rhoinv = hRep.rho_inverse_[node];
		dPhi = multStateAB(rhoinv, dPhi);
	}

	/// 4.) Multiply with propagation phase (for real- or imaginary time propagation)
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

////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////

void symLayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase) {

	dPhi = TreeFunctions::symApply(Phi, hRep.hMatSets_, H, node);
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
