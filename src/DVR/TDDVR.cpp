//
// Created by Roman Ellerbrock on 3/9/20.
//

#include "TDDVR.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Util/WeightedSimultaneousDiagonalization.h"

vector<Matrixcd> getXs(const vector<SparseMatrixTreecd>& Xs, const Node& node) {
	vector<Matrixcd> xs;
	for (const auto& xtree : Xs) {
		if (xtree.Active(node)) {
			xs.push_back(xtree[node]);
		}
	}
	return xs;
}

void setGrids(TreeGrids& gridtrees, vector<Vectord> grids, const Node& node) {
	size_t k = 0;
	for (auto& grid : gridtrees) {
		if (grid.Active(node)) {
			grid[node] = grids[k];
			k++;
		}
	}
}

Matrixcd Regularize(Matrixcd w) {
	double norm = abs(w.Trace());
	w = (1. / norm) * w;
	double eps = 1e-6;
	for (size_t i = 0; i < w.Dim1(); ++i) {
		w(i, i) += eps * exp(-w(i, i) / eps);
	}
	// Normalize the weight matrix
	norm = abs(w.Trace());
	w = (1. / norm) * w;
	return w;
}

void LayerGrid(TreeGrids& grids, Matrixcd& trafo, const vector<SparseMatrixTreecd>& Xs,
	const Matrixcd* w_ptr, const Node& node) {
	assert(grids.size() == Xs.size());
	auto xs = getXs(Xs, node);
	Matrixcd w;
	if (w_ptr != nullptr) {
		w = *w_ptr;
		w = Regularize(w);
	} else {
		w = IdentityMatrix<complex<double>>(trafo.Dim1());
	}

	assert(!xs.empty());
	auto diags = WeightedSimultaneousDiagonalization::Calculate(xs, w, 1e-12);
	setGrids(grids, diags.second, node);
	trafo = diags.first;
}

void TDDVR::Calculate(const Wavefunction& Psi, const Tree& tree) {
	/// Calculate density matrix
	TreeFunctions::Contraction(rho_, Psi, tree, true);

	/// Calculate X-Matrices
	Xs_.Update(Psi, tree);

	/// Build standard grid
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			LayerGrid(grids_, trafo_[node], Xs_.xmats_, &rho_[node], node);
		}
	}

	/// Build hole grid
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			LayerGrid(hole_grids_, hole_trafo_[node], Xs_.xholes_, nullptr, node);
		}
	}
}

