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

class SymXMatrixTrees {
public:

	SymXMatrixTrees(const Tree& tree)
	: xops_(Xsop(tree)) {
		xmat_.clear();
		for (const auto& x : xops_) {
			SparseMatrixTreecd x1(x, tree);
			SparseMatrixTreecd x2(x, tree);
			SparseMatrixTreePaircd y({x1, x2});
			xmat_.emplace_back(y);
		}
	}

	~SymXMatrixTrees() = default;

	void update(const SymTensorTree& Psi, const Tree& tree) {
		TreeFunctions::symRepresent(xmat_, Psi, Psi, xops_, tree);
	}

	SOPcd xops_;

	SymMatrixTrees xmat_;
};


#endif //SYMXMATRIXTREES_H
