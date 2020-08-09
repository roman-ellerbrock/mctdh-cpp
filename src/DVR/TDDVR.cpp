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

vector<double> CalculateShift(const vector<Matrixcd>& xs, const Matrixcd& w) {
	vector<double> shifts;
	for (const auto& x : xs) {
		double s = abs(x.Trace());
		shifts.push_back(s);
	}
	return shifts;
}

void shift(vector<Matrixcd>& xs, const vector<double>& shift) {
	assert(xs.size() == shift.size());
	for (size_t i = 0; i < xs.size(); ++i) {
		double s = shift[i];
		Matrixcd& x = xs[i];
		for (size_t j = 0; j < x.Dim1(); ++j) {
			x(j, j) -= s;
		}
	}
}

void shift(vector<Vectord>& xs, const vector<double>& shift) {
	assert(xs.size() == shift.size());
	for (size_t i = 0; i < xs.size(); ++i) {
		double s = shift[i];
		Vectord& x = xs[i];
		for (size_t j = 0; j < x.Dim(); ++j) {
			x(j) += s;
		}
	}
}

void LayerGrid(TreeGrids& grids, Matrixcd& trafo,
	const vector<SparseMatrixTreecd>& Xs,
	const Matrixcd *w_ptr, const Node& node) {
	assert(grids.size() == Xs.size());
	auto xs = getXs(Xs, node);
	Matrixcd w;
	if (w_ptr != nullptr) {
		w = *w_ptr;
		w = Regularize(w, 1e-8);
	} else {
		w = IdentityMatrix<complex<double>>(trafo.Dim1());
	}

	assert(!xs.empty());

	auto shifts = CalculateShift(xs, w);

	shift(xs, shifts);
/*	node.info();
	for (auto s : shifts) {
		cout << s << "\t";
	}
	cout << endl;*/
	auto diags = WeightedSimultaneousDiagonalization::Calculate(xs, w, 1e-12);
	shift(diags.second, shifts);
	setGrids(grids, diags.second, node);
	trafo = diags.first;
	trafo = trafo.Adjoint();
}

void UpdateGrids(TreeGrids& grids, MatrixTreecd& trafo, const vector<SparseMatrixTreecd>& Xs,
	const MatrixTreecd *rho_ptr, const Tree& tree) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			if ((rho_ptr == nullptr) || node.isBottomlayer()) {
//			if ((rho_ptr == nullptr) ) {
				LayerGrid(grids, trafo[node], Xs, nullptr, node);
			} else {
				LayerGrid(grids, trafo[node], Xs, &rho_ptr->operator[](node), node);
			}
		}
	}
}

void TDDVR::Update(const Wavefunction& Psi, const Tree& tree) {
	/// Calculate density matrix
	TreeFunctions::Contraction(rho_, Psi, tree, true);

	/// Calculate X-Matrices
	Xs_.Update(Psi, tree);

	/// Build standard grid
	UpdateGrids(grids_, trafo_, Xs_.mats_, &rho_, tree);

	/// Build hole grid
	UpdateGrids(hole_grids_, hole_trafo_, Xs_.holes_, &rho_, tree);
}

void TDDVR::GridTransformationLocal(Tensorcd& Phi, const Node& node, bool inverse) const {

	/// Transform underlying A-coefficient
	for (size_t k = 0; k < node.nChildren(); ++k) {
		const Node& child = node.child(k);
		if (!inverse) {
			Phi = MatrixTensor(trafo_[child], Phi, k);
		} else {
			Phi = multATB(trafo_[child], Phi, k);
		}
	}

	/// Transform state
	if (!inverse) {
		Phi = multStateArTB(hole_trafo_[node], Phi);
	} else {
		auto m = hole_trafo_[node].Adjoint();
		Phi = multStateArTB(m, Phi);
	}
}

void TDDVR::GridTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node : tree) {
		GridTransformationLocal(Psi[node], node, inverse);
	}
}

void TDDVR::NodeTransformation(Tensorcd& Phi, const Node& node, bool inverse) const {

	/// Transform underlying A-coefficient
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			if (!inverse) {
				Phi = MatrixTensor(trafo_[child], Phi, k);
			} else {
				Phi = multATB(trafo_[child], Phi, k);
			}
		}
	}

	/// Transform state
	if (!node.isToplayer()) {
		if (!inverse) {
			Phi = MatrixTensor(hole_trafo_[node], Phi, node.nChildren());
		} else {
			Phi = multATB(hole_trafo_[node], Phi, node.nChildren());
		}
	}
}

void TDDVR::EdgeTransformation(Matrixcd& B_inv, const Edge& edge, bool inverse) const {
	if (!inverse) {
		B_inv = trafo_[edge].Transpose().Adjoint() * B_inv;
		B_inv = B_inv * hole_trafo_[edge].Adjoint();
	} else {
		B_inv = trafo_[edge].Transpose() * B_inv;
		B_inv = B_inv * hole_trafo_[edge];
	}
}

void TDDVR::NodeTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node : tree) {
		NodeTransformation(Psi[node], node, inverse);
	}
}

void TDDVR::EdgeTransformation(MatrixTreecd& B_inv, const Tree& tree, bool inverse) const {
	for (const Edge& edge : tree.Edges()) {
		EdgeTransformation(B_inv[edge], edge, inverse);
	}
}

void TDDVR::GridTransformation(ExplicitEdgeWavefunction& Psi, const Tree& tree, bool inverse) const {
	NodeTransformation(Psi.nodes(), tree, inverse);
	EdgeTransformation(Psi.edges(), tree, inverse);
}

void TDDVR::print(const Tree& tree) const {
	cout << "TDDVR: " << endl;
	cout << "Grids:" << endl;
	for (const Node& node : tree) {
		if (!node.isToplayer() && !node.isBottomlayer()) {
			size_t dim = trafo_[node].Dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const SparseVectorTreed& grid : grids_) {
					if (grid.Active(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
	}
	cout << "Hole grids:" << endl;
	for (const Node& node : tree) {
		if (!node.isToplayer() && !node.isBottomlayer()) {
			size_t dim = trafo_[node].Dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const SparseVectorTreed& grid : hole_grids_) {
					if (grid.Active(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
	}
}




