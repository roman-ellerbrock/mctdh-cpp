//
// Created by Roman Ellerbrock on 1/4/21.
//

#ifndef SYMXMATRIXTREES_H
#define SYMXMATRIXTREES_H
#include "MatrixTensorTree.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "MatrixTensorTreeFunctions.h"
#include "XMatrixTrees.h"
#include "TreeClasses/SymTensorTree.h"

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
			auto x2 = SparseMatrixTreecd(M, tree, false, false);
//			auto x2 = SparseMatrixTreecd(M, tree, false, true);
			SparseMatrixTreePaircd y({x1, x2});
			emplace_back(y);
		}
	}

	~SymXMatrixTrees() = default;

	void Update(const Wavefunction& Psi, const Tree& tree) {
		cout << "Casting to sTTN:\n";
		SymTensorTree psi(Psi, tree);
		psi.up_.print(tree);

		cout << "representing:\n";
		for (size_t l = 0; l < xops_.size(); ++l) {
			TreeFunctions::represent(xmat_[l].first, xops_[l], psi.up_, psi.up_, tree);
		}

//		TreeFunctions::symRepresent(xmat_, psi, psi, xops_, tree);

		size_t l = 0;
		for (const auto& xmat : xmat_) {
			cout << "coordinate l = " << l <<endl;
			cout << "mat:\n";
			xmat.first.print();
			cout << "hole:\n";
			xmat.second.print();
		}
	}

	void Optimize(MatrixTensorTree& Psi, const Tree& tree);

	MatrixTensorTree OptimizeUp(const MatrixTensorTree& Psi,
		const Tree& tree, const Tree& tree_small) const;

	Tensorcd OptimizeUp(const Tensorcd& Phi, const Matrixcd& rho,
		const Node& node, const Node& node_small) const;

	SOPcd xops_;

	SymMatrixTrees xmat_;
};


#endif //SYMXMATRIXTREES_H
