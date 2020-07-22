//
// Created by Roman Ellerbrock on 3/11/20.
//
#include "DVR/XMatrixTrees.h"
#include <random>

Wavefunction Regularize(Wavefunction Psi, const Tree& tree, double eps) {
	const Node& top = tree.TopNode();
	const Tensorcd& A = Psi[top];
	const TensorShape& shape = A.shape();
	normal_distribution<double> dist;
	mt19937 gen(23949);
	for (size_t i = 0; i < shape.totalDimension(); ++i) {
		A(i) += eps * dist(gen);
	}
	return Psi;
}

SOPcd Xsop(const Tree& tree) {
	LeafFuncd x = &LeafInterface::applyX;
	LeafFuncd I = &LeafInterface::Identity;
	SOPcd xops;
	for (size_t l = 0; l < tree.nLeaves(); ++l) {
		const Leaf& leaf = tree.GetLeaf(l);
		size_t mode = leaf.Mode();
		MLOcd M(x, mode);
		for (size_t i = 0; i < tree.nLeaves(); ++i) {
			M.push_back(I, i);
		}
		xops.push_back(M, 1.);
	}
	return xops;
}

Tensorcd XMatrixTrees::Optimize(const Tensorcd& Phi, const Matrixcd& rho,
	const Node& node, const Node& node_small) const {

	size_t n_occ = node_small.shape().lastDimension();

	auto X = BuildX(Phi, rho, node);
	X = UnProject(n_occ, X, Phi);
	auto xspec = Diagonalize(X);
	auto oPhi = Occupy(Phi, xspec.first, n_occ, node);

	return oPhi;
}

Wavefunction XMatrixTrees::Optimize(Wavefunction Psi,
	const MatrixTreecd& rho, const Tree& tree, const Tree& tree_small) {

	Wavefunction Chi(Psi);
	for (const Node& node : tree) {
		Chi[node] = Optimize(Psi[node], rho[node],
			node, tree_small.GetNode(node.Address()));
	}
	return Chi;
}

Matrixcd XMatrixTrees::BuildX(const Tensorcd& Phi, const Matrixcd& rho,
	const Node& node) const {

	const Leaf& leaf = node.getLeaf();
	const LeafInterface& grid = leaf.PrimitiveGrid();
	auto w = 1. / abs(rho.Trace()) * rho;
//	w = Regularize(w, 1e-5);
	Tensorcd xPhi(Phi.shape());
	grid.applyX(xPhi, Phi);
	auto wxPhi = MatrixTensor(w, xPhi, node.parentIdx());
	return Tensor_Extension::OuterProduct(wxPhi, Phi);
}

Matrixcd UnProject(size_t n_occupied, const Matrixcd& X,
	const Tensorcd& Phi) {
	/**
	 * Calculate (1 - P) A (1 - P) from A
	 */
	size_t dimpart = Phi.shape().lastBefore();
	assert (Phi.shape().lastDimension() >= n_occupied);

	// Build unit
	Matrixcd Pu(dimpart, dimpart);
	for (size_t i = 0; i < dimpart; ++i)
		Pu(i, i) = 1;

	// Build P (it doesnt matter if you use ntesor or small_ntensor
	for (size_t i = 0; i < dimpart; ++i) {
		for (size_t j = 0; j < dimpart; ++j) {
			for (size_t n = 0; n < n_occupied; ++n) {
				Pu(j, i) -= Phi(j, n) * conj(Phi(i, n));
			}
		}
	}
	// (1 - P) A (1 - P)
	Matrixcd B = X * Pu;
	return Pu * B;
}

Tensorcd Occupy(const Tensorcd& Phi, const Matrixcd& trafo,
	size_t n_occupied, const Node& node) {

	const TensorShape& shape = node.shape();
	size_t dimpart = shape.lastBefore();
	size_t ntensor = shape.lastDimension();
	assert(ntensor >= n_occupied);

	// Copy the highest value eigenvectors to a new Phi
	Tensorcd oPhi(Phi);

	// (The following ints have to be integers, not size_t)
	int last = (int) dimpart - 1;
	for (int n = 0; n < (ntensor - n_occupied); ++n) {
		int mat_idx = last - n;
		int t_idx = (int) n_occupied + n;
		for (size_t i = 0; i < dimpart; i++) {
			oPhi(i, (size_t) t_idx) = trafo(i, (size_t) mat_idx);
		}
	}

	GramSchmidt(oPhi);
	return oPhi;
}

