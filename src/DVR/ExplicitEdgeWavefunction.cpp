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
	edges_ = sqrt(rho, tree);
	for (auto& e : edges_) { e = e.Transpose(); }

	/// Build node representation (A^\tilde's)
	for (const Edge& e : tree.Edges()) {
		const Node& node = e.down();
		const Matrixcd& B = edges_[e];
		Tensorcd& A = nodes_[node];
		/// Basically like multiplying with sqrt(rho)'s
		A = TensorMatrix(A, B, node.nChildren());
	}

	B_inv_ = inverse(edges(), tree);
}

TensorTreecd ExplicitEdgeWavefunction::TopDownNormalized(const Tree& tree) const {
	/// Contraction-normalized representation
	/// Build A^{(p\circ k) p}
	/// Note: Tensors get moved one layer down!
	TensorTreecd Psi(nodes());
	for (const Edge& e : tree.Edges()) {
		const Node& node = e.down();
		const Node& parent = node.parent();
		Psi[node] = TensorMatrix(Psi[parent], B_inv_[e], node.childIdx());
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
		Psi[node] = TensorMatrix(Psi[node], B_inv_[e], node.nChildren());
	}

	return Psi;
}
