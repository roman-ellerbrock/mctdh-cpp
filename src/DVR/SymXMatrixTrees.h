//
// Created by Roman Ellerbrock on 1/4/21.
//

#ifndef SYMXMATRIXTREES_H
#define SYMXMATRIXTREES_H
#include "MatrixTensorTree.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "MatrixTensorTreeFunctions.h"
#include "XMatrixTrees.h"

class SymXMatrixTrees : public SparseMatrixTreePairscd {
public:

	SymXMatrixTrees(const Tree& tree)
	: xops_(Xsop(tree)) {
		LeafFuncd x = &LeafInterface::applyX;

		clear();
		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.getLeaf(l);
			size_t mode = leaf.mode();
			MLOcd M(x, mode);

			auto x1 = SparseMatrixTreecd(M, tree);
			auto x2 = SparseMatrixTreecd(M, tree, false, true);
			SparseMatrixTreePaircd y({x1, x2});
			emplace_back(y);
		}
	}

	~SymXMatrixTrees() = default;

	void Update();

	void Optimize(MatrixTensorTree& Psi, const Tree& tree);

	MatrixTensorTree OptimizeUp(const MatrixTensorTree& Psi,
		const Tree& tree, const Tree& tree_small) const;

	Tensorcd OptimizeUp(const Tensorcd& Phi, const Matrixcd& rho,
		const Node& node, const Node& node_small) const;

	SOPcd xops_;
};


#endif //SYMXMATRIXTREES_H
