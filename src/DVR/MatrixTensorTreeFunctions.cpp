//
// Created by Roman Ellerbrock on 1/3/21.
//

#include "MatrixTensorTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"

namespace TreeFunctions {

	void Represent(SparseMatrixTreecd& mat, const MatrixTensorTree& Psi, const MLOcd& M,
		const Tree& tree) {

		const TensorTreecd& Psi_up = Psi.BottomUpNormalized(tree);
		TreeFunctions::Represent(mat, M, Psi_up, tree);
	}

	void ContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Phi,
		const SparseMatrixTreecd& mat, const Node& hchild) {
		assert(!hchild.isToplayer());
		const Node& parent = hchild.parent();

		Tensorcd hPhi = TreeFunctions::ApplyHole(hole, Phi, hchild);

		if (!parent.isToplayer()) {
			hPhi = TensorMatrix(hPhi, hole[parent], parent.childIdx());
		}
		hole[hchild] = Contraction(Phi, hPhi, hchild.childIdx());
	}

	void Contraction(SparseMatrixTreecd& hole, const MatrixTensorTree& Psi,
		const SparseMatrixTreecd& mat, const SparseTree& marker, const Tree& tree) {

		const TensorTreecd& Psi_down = Psi.TopDownNormalized(tree);
		int sub_topnode = marker.size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = marker.MCTDHNode(n);
			if (!node.isToplayer()) {
				ContractionLocal(hole, Psi_down[node], mat, node);
			}
		}

	}

	typedef pair<SparseMatrixTreecd, SparseMatrixTreecd> SparseMatrixTreePaircd;
	void Represent(SparseMatrixTreePaircd& mats,
		const MatrixTensorTree& Psi, const MLOcd& M,
		const SparseTree& stree, const Tree& tree) {

		Represent(mats.first, Psi, M, tree);
		Contraction(mats.second, Psi, mats.first, stree, tree);

	}
}