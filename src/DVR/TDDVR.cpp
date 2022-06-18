//
// Created by Roman Ellerbrock on 3/9/20.
//

#include <Core/TensorBLAS.h>
#include "TDDVR.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Util/WeightedSimultaneousDiagonalization.h"
#include "Util/SimultaneousDiagonalization.h"

vector<Matrixcd> getXs(const vector<SparseMatrixTreecd>& Xs, const Node& node) {
	vector<Matrixcd> xs;
	for (const auto& xtree : Xs) {
		if (xtree.isActive(node)) {
			xs.push_back(xtree[node]);
		}
	}
	return xs;
}

void setGrids(TreeGrids& gridtrees, vector<Vectord> grids, const Node& node) {
	size_t k = 0;
	for (auto& grid : gridtrees) {
		if (grid.isActive(node)) {
			grid[node] = grids[k];
			k++;
		}
	}
}

vector<double> calculateShift(const vector<Matrixcd>& xs, const Matrixcd& w) {
	vector<double> shifts;
	for (const auto& x : xs) {
		double s = abs(x.trace());
		shifts.push_back(s);
	}
	return shifts;
}

void shift(vector<Matrixcd>& xs, const vector<double>& shift) {
	assert(xs.size() == shift.size());
	for (size_t i = 0; i < xs.size(); ++i) {
		double s = shift[i];
		Matrixcd& x = xs[i];
		for (size_t j = 0; j < x.dim1(); ++j) {
			x(j, j) -= s;
		}
	}
}

void shiftBack(vector<Vectord>& xs, const vector<double>& shift) {
	assert(xs.size() == shift.size());
	for (size_t i = 0; i < xs.size(); ++i) {
		double s = shift[i];
		Vectord& x = xs[i];
		for (size_t j = 0; j < x.dim(); ++j) {
			x(j) += s;
		}
	}
}

void LayerGrid(TreeGrids& grids, Matrixcd& trafo,
	const vector<SparseMatrixTreecd>& Xs,
	const Matrixcd *w_ptr, const Node& node) {
	assert(grids.size() == Xs.size());
	auto xs = getXs(Xs, node);
	assert(!xs.empty());

	if (xs.size() == 1) {
		SpectralDecompositioncd diags = diagonalize(xs.front());
		trafo = diags.first;
		trafo = trafo.adjoint();
		setGrids(grids, {diags.second}, node);
	} else {
		Matrixcd w;
		if (w_ptr != nullptr) {
			w = *w_ptr;
			w = regularize(w, 1e-9);
		} else {
			w = identityMatrix<complex<double>>(trafo.dim1());
		}

		auto shifts = calculateShift(xs, w);
		shift(xs, shifts);

		auto diags = WeightedSimultaneousDiagonalization::calculate(xs, w, 1e-5);

		shiftBack(diags.second, shifts);

		setGrids(grids, diags.second, node);
		trafo = diags.first;
		trafo = trafo.adjoint();
	}
}

void UpdateGrids(TreeGrids& grids, MatrixTreecd& trafo, const vector<SparseMatrixTreecd>& Xs,
	const MatrixTreecd *rho_ptr, const Tree& tree) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
//			if ((rho_ptr == nullptr) || node.isBottomlayer()) {
			if (rho_ptr == nullptr) {
				LayerGrid(grids, trafo[node], Xs, nullptr, node);
			} else {
				LayerGrid(grids, trafo[node], Xs, &rho_ptr->operator[](node), node);
			}
		}
	}
}

/*
#include "Util/RandomMatrices.h"

void wavelet(const Tree& tree) {
	for (size_t l = 0; l < tree.nLeaves(); ++l) {
		const Leaf& leaf = tree.getLeaf(l);
		LeafFuncd xl = &LeafInterface::applyX;
		LeafFuncd pl = &LeafInterface::applyP;
		LeafFunctioncd xf(xl);
		LeafFunctioncd pf(pl);

		vector<Matrixcd> xps;
		xps.push_back(toMatrix(xf, leaf));
		xps.push_back(toMatrix(pf, leaf));
		Matrixcd& x = xps.front();
		Matrixcd& p = xps.back();

		Matrixcd trafo(xps.front());
		trafo.zero();
		for (size_t i = 0; i < trafo.dim1(); ++i) {
			trafo(i, i) = 1.;
		}
		size_t dim = trafo.dim1();


		/// random init
		mt19937 gen(time(nullptr));
		auto U = RandomMatrices::gue(trafo.dim1(), gen);
//		x = U * x * U.adjoint();
//		p = U * p * U.adjoint();


		cout << "Leaf: " << l << endl;
		double normx = x.frobeniusNorm();
		double normp = p.frobeniusNorm();
		x /= normx;
		p /= normp;
		cout << "norm(x) = " << x.frobeniusNorm() << endl;
		cout << "norm(p) = " << p.frobeniusNorm() << endl;
		auto comm = x * p - p * x;
		cout << "[x, p] = " << comm.frobeniusNorm() << endl;

		auto xx = diagonalize(x);
		auto px = diagonalize(p);

		SimultaneousDiagonalization sd;
		sd.initialization(xps, 1e-10);
		sd.calculate(xps, trafo);

		Vectord xvec(dim);
		Vectord pvec(dim);
		for (size_t i = 0; i < dim; ++i) {
			xvec(i) = real(x(i, i));
			pvec(i) = real(p(i, i));
		}

		cout << "x:\n";
		xvec.print();
		cout << endl;
		cout << "p:\n";
		pvec.print();
		cout << endl;

		cout << "x diag:\n";
		xx.second.print();
		cout << endl;

		cout << "p diag:\n";
		px.second.print();
		cout << endl;
		cout << "Example: " << xvec(20) << endl;
		for (size_t i = 0; i < dim; ++i) {
			cout << pow(abs(trafo(i, 40)), 2) << " " << pow(abs(xx.first(i, 40)), 2) << endl;
		}
		getchar();
	}
}
*/

void TDDVR::Update(const Wavefunction& Psi, const Tree& tree) {
	/// Calculate density matrix
	TreeFunctions::contraction(rho_, Psi, tree, true);

	/// Calculate X-Matrices
	Xs_.Update(Psi, tree);

	/// Build standard grid
	UpdateGrids(grids_, trafo_, Xs_.mats_, &rho_, tree);

	/// Build hole grid
	UpdateGrids(hole_grids_, hole_trafo_, Xs_.holes_, &rho_, tree);
}

void TDDVR::NodeTransformation(Tensorcd& Phi, const Node& node, bool inverse) const {

	/// @TODO: make in/out independent
	/// Transform underlying A-coefficient
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
//			Tensorcd Xi(Phi.shape(), &(mem_.work1_[0]), false, false);
			if (!inverse) {
				Phi = matrixTensorBLAS(trafo_[child], Phi, k);
//				matrixTensorBLAS(Xi, mem_.work2_, trafo_[child], Phi, k, true);
			} else {
				Phi = matrixTensorBLAS(trafo_[child].adjoint(), Phi, k);
//				matrixTensorBLAS(Xi, mem_.work2_, trafo_[child].adjoint(), Phi, k, true);
			}
//			Phi = Xi;
		}
	}

	/// Transform state
	if (!node.isToplayer()) {
//		Tensorcd Xi(Phi.shape(), &(mem_.work1_[0]), false, false);
		if (!inverse) {
			Phi = matrixTensorBLAS(hole_trafo_[node], Phi, node.nChildren());
//			Xi.shape().print();
//			matrixTensorBLAS(Xi, mem_.work2_, hole_trafo_[node], Phi, node.parentIdx(), true);
		} else {
			Phi = matrixTensorBLAS(hole_trafo_[node].adjoint(), Phi, node.nChildren());
//			matrixTensorBLAS(Xi, mem_.work2_, hole_trafo_[node].adjoint(), Phi, node.parentIdx(), true);
		}
//		Phi = Xi;
	}
}

void TDDVR::EdgeTransformation(Matrixcd& B_inv, const Edge& edge, bool inverse) const {
	if (!inverse) {
		B_inv = trafo_[edge].transpose().adjoint() * B_inv;
		B_inv = B_inv * hole_trafo_[edge].adjoint();
	} else {
		B_inv = trafo_[edge].transpose() * B_inv;
		B_inv = B_inv * hole_trafo_[edge];
	}
}

void TDDVR::NodeTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node : tree) {
		NodeTransformation(Psi[node], node, inverse);
	}
}

void TDDVR::EdgeTransformation(MatrixTreecd& B_inv, const Tree& tree, bool inverse) const {
	for (const Edge& edge : tree.edges()) {
		EdgeTransformation(B_inv[edge], edge, inverse);
	}
}

void TDDVR::GridTransformation(MatrixTensorTree& Psi, const Tree& tree, bool inverse) const {
	NodeTransformation(Psi.nodes(), tree, inverse);
	EdgeTransformation(Psi.edges(), tree, inverse);
}

void TDDVR::print(const Tree& tree) const {
	cout << "TDDVR: " << endl;
	cout << "Grids:" << endl;
	for (const Node& node : tree) {
		if (!node.isToplayer() && !node.isBottomlayer()) {
			size_t dim = trafo_[node].dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const SparseVectorTreed& grid : grids_) {
					if (grid.isActive(node)) {
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
			size_t dim = trafo_[node].dim1();
			node.info();
			for (size_t i = 0; i < dim; ++i) {
				for (const SparseVectorTreed& grid : hole_grids_) {
					if (grid.isActive(node)) {
						const Vectord& g = grid[node];
						cout << g(i) << "\t";
					}
				}
				cout << endl;
			}
		}
	}
}


/*void TDDVR::GridTransformationLocal(Tensorcd& Phi, const Node& node, bool inverse) const {

	/// Transform underlying A-coefficient
	for (size_t k = 0; k < node.nChildren(); ++k) {
		const Node& child = node.child(k);
		if (!inverse) {
			Phi = matrixTensorBLAS(trafo_[child], Phi, k);
		} else {
			Phi = matrixTensorBLAS(trafo_[child].adjoint(), Phi, k);
		}
	}

	/// Transform state
	if (!inverse) {
//		Phi = tensorMatrix(Phi, hole_trafo_[node], node.parentIdx());
		Phi = matrixTensorBLAS(hole_trafo_[node].transpose(), Phi, node.parentIdx());
	} else {
//		Phi = tensorMatrix(Phi, hole_trafo_[node].adjoint(), node.parentIdx());
		Phi = matrixTensorBLAS(hole_trafo_[node].adjoint().transpose(), Phi, node.parentIdx());
	}
}*/

/*void TDDVR::GridTransformation(Wavefunction& Psi, const Tree& tree, bool inverse) const {
	for (const Node& node : tree) {
		GridTransformationLocal(Psi[node], node, inverse);
	}
}*/



