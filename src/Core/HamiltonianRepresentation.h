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
#include "DVR/CDVR.h"
#include "DVR/MatrixTensorTreeFunctions.h"

#include "TreeShape/LeafTypes/HO_Basis.h"
#include "TreeShape/LeafTypes/FFTGrid.h"
#include "TreeShape/LeafTypes/LegendrePolynomials.h"

class HamiltonianRepresentation {
public:
	HamiltonianRepresentation(const Hamiltonian& H, const Tree& tree,
		const Tree& cdvrtree)
		: rho_(tree), rho_decomposition_(tree),
		rho_inverse_(tree), cdvr_(tree), mem_(tree) {
		;
		hMats_.clear();
		hContractions_.clear();
		for (const auto& M : H) {
			hMats_.emplace_back(SparseMatrixTreecd(M, tree, true));
			hContractions_.emplace_back(SparseMatrixTreecd(M, tree));
		}

		hMatSets_.clear();
		for (const auto& M : H) {
			auto x1 = SparseMatrixTreecd(M, tree, false);
			auto x2 = SparseMatrixTreecd(M, tree);
			SparseMatrixTreePaircd y({x1, x2});
			hMatSets_.emplace_back(y);
		}
	}

	HamiltonianRepresentation(const Hamiltonian& H, const Tree& tree)
		: HamiltonianRepresentation(H, tree, tree) {
	}

	~HamiltonianRepresentation() = default;

	void build(const Hamiltonian& H, const Wavefunction& Psi, const Tree& tree, double time) {

		/// Calculate density matrix tree
		TreeFunctions::contraction(rho_, Psi, tree, true);

		/// Density matrix tree decomposition
		rho_decomposition_.calculate(rho_, tree);

		/// Calculate inverse density matrix tree
		rho_inverse_ = rho_decomposition_.invert(tree, 1e-10);
//		rho_inverse_ = inverse(rho_, tree);

		/// Calculate h-matrix trees
		TreeFunctions::represent(hMats_, H, Psi, Psi, tree);

		/// Calculate h-matrix tree contractions
		TreeFunctions::contraction(hContractions_, hMats_, Psi, Psi, tree);

		/// Calculate CDVR
		if (H.hasV) { cdvr_.Update(Psi, H.V_, tree); }
//		if (H.hasV) { cdvr_.Update2(Psi, H.V_, tree); }
	}

	void symbuild(const Hamiltonian& H, MatrixTensorTree Psi,
		const Tree& tree) {

		Psi.buildFromWeighted(tree);

		TreeFunctions::Represent(hMatSets_, Psi, H, tree);
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

	int counter_ = 0;

	MatrixTreecd rho_;
	SpectralDecompositionTreecd rho_decomposition_;
	MatrixTreecd rho_inverse_;
	SparseMatrixTreescd hMats_;
	SparseMatrixTreescd hContractions_;

	SparseMatrixTreePairscd hMatSets_;

	WorkMemorycd mem_;

	CDVR cdvr_;
};

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep, const Node& node);

Matrixcd Expectation(const HamiltonianRepresentation& hRep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree);

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase = 1.);

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase = 1.);

void symDerivative(MatrixTensorTree& dPsi, HamiltonianRepresentation& hRep,
	double time, const MatrixTensorTree& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase = 1.);

#endif //HAMILTONIANREPRESENTATION_H
