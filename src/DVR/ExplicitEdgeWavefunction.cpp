//
// Created by Roman Ellerbrock on 4/3/20.
//

#include "ExplicitEdgeWavefunction.h"
#include "TreeClasses/SpectralDecompositionTree.h"

ExplicitEdgeWavefunction::ExplicitEdgeWavefunction(const Wavefunction& Psi, const Tree& tree, bool orthogonal) {
	/// Note: Requires orthogonal wavefunction representation (typically given)
	assert(orthogonal);
	Wavefunction& nodes_ = first;
	MatrixTreecd& edges_ = second;
	nodes_ = Psi;

	/// Get edge matrices (B's)
	MatrixTreecd rho = TreeFunctions::Contraction(nodes_, tree, true);
	auto B = sqrt(rho, tree);

	/// Build node representation (A^\tilde's)
	for (const Edge& e : tree.Edges()) {
		const Node& node = e.down();
		Tensorcd& A = nodes_[node];
		/// Basically like multiplying with sqrt(rho)'s
		A = MatrixTensor(B[e], A, node.nChildren());
	}

	edges_= inverse(B, tree);
}

TensorTreecd ExplicitEdgeWavefunction::TopDownNormalized(const Tree& tree) const {
	/// Contraction-normalized representation
	/// Build A^{(p\circ k) p}
	/// Note: Tensors get moved one layer down!
	TensorTreecd Psi(nodes());
	for (const Edge& e : tree.Edges()) {
		const Node& node = e.down();
		const Node& parent = node.parent();
		Psi[node] = MatrixTensor(edges()[e].Transpose(), Psi[parent], node.childIdx());
	}

	return Psi;
}

TensorTreecd ExplicitEdgeWavefunction::BottomUpNormalized(const Tree& tree) const {
	/// Dot-Product normalized reperesentation
	/// Build A^{p\circ k (p)}
	/// This is the conventional wavefunction representation.
	TensorTreecd Psi(nodes());

	for (const Edge& e : tree.Edges()) {
		const Node& node = e.down();
		Psi[node] = MatrixTensor(edges()[e], Psi[node], node.nChildren());
	}

	return Psi;
}

bool IsWorking_bottomup(const ExplicitEdgeWavefunction& Psi, const Tree& tree, double eps) {
	auto bottomup = Psi.BottomUpNormalized(tree);
	for (const Node& node : tree) {
		auto x = Contraction(bottomup[node], bottomup[node], node.nChildren());
		auto r = Residual(x, IdentityMatrixcd(x.Dim1()));
		if (r > eps) { cerr << "bottom-up normalization failed.\n"; return false; }
	}
	return true;
}

bool IsWorking_topdown(const ExplicitEdgeWavefunction& Psi, const Tree& tree, double eps) {
	auto topdown = Psi.TopDownNormalized(tree);
	for (const Edge& e : tree.Edges()) {
		const Node& node = e.down();
		auto x = Contraction(topdown[node], topdown[node], node.childIdx());
		auto r = Residual(x, IdentityMatrixcd(x.Dim1()));
		if (r > eps) { cerr << "top-down normalization failed.\n"; return false; }
	}
	return true;
}

bool IsWorking(const ExplicitEdgeWavefunction& Psi, const Tree& tree, double eps) {
	return (IsWorking_bottomup(Psi, tree, eps) && IsWorking_topdown(Psi, tree, eps));
}
