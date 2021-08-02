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

template<typename T>
class tHamiltonianRepresentation {
public:
	tHamiltonianRepresentation(const SOP<T>& H, const Tree& tree)
		: rho_(tree), rho_decomposition_(tree), rho_inverse_(tree) {
		hMats_.clear();
		hContractions_.clear();
		for (const auto& M : H) {
			hMats_.emplace_back(SparseMatrixTree<T>(M, tree, true));
			hContractions_.emplace_back(SparseMatrixTree<T>(M, tree));
		}
	}

	~tHamiltonianRepresentation() = default;

	void build(const SOP<T>& H, const TensorTree<T>& Psi, const Tree& tree) {
		/// Calculate density matrix tree
		TreeFunctions::contraction(rho_, Psi, tree, true);

		/// Density matrix tree decomposition
		rho_decomposition_.calculate(rho_, tree);

		/// Calculate inverse density matrix tree
		rho_inverse_ = rho_decomposition_.invert(tree);
//		rho_inverse_ = inverse(rho_, tree);

		/// Calculate h-matrix trees
		TreeFunctions::represent(hMats_, H, Psi, Psi, tree);

		/// Calculate h-matrix tree contractions
		TreeFunctions::contraction(hContractions_, hMats_, Psi, Psi, tree);
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

	MatrixTree<T> rho_;
	SpectralDecompositionTree<T> rho_decomposition_;
	MatrixTree<T> rho_inverse_;
	SparseMatrixTrees<T> hMats_;
	SparseMatrixTrees<T> hContractions_;
};

template<typename T>
Tensorcd Apply(const SOP<T>& H, const Tensor<T>& Phi,
	const tHamiltonianRepresentation<T>& hRep, const Node& node);

template<typename T>
Matrixcd Expectation(const tHamiltonianRepresentation<T>& hRep,
	const TensorTree<T>& Psi, const SOP<T>& H, const Tree& tree);

template<typename T>
void LayerDerivative(Tensor<T>& dPhi, double time, const Tensor<T>& Phi,
	const SOP<T>& H, const tHamiltonianRepresentation<T>& hRep,
	const Node& node, T propagation_phase = 1.);

template<typename T>
void Derivative(TensorTree<T>& dPsi, tHamiltonianRepresentation<T>& hRep,
	double time, const TensorTree<T>& Psi, const SOP<T>& H,
	const Tree& tree, T propagation_phase = 1.);

#endif //HAMILTONIANREPRESENTATION_H
