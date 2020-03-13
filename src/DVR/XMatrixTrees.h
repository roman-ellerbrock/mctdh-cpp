//
// Created by Roman Ellerbrock on 3/11/20.
//

#ifndef XMATRIXTREES_H
#define XMATRIXTREES_H
#include "Core/Wavefunction.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/SOPMatrixTrees.h"

SOPcd Xsop(const Tree& tree);

class XMatrixTrees {
public:
	explicit XMatrixTrees(const Tree& tree)
		: xops_(Xsop(tree)) {
		LeafFuncd x = &LeafInterface::applyX;

		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			size_t mode = leaf.Mode();
			MLOcd M(x, mode);

			xmats_.emplace_back(SparseMatrixTreecd(M, tree));
			xholes_.emplace_back(SparseMatrixTreecd(M, tree, true, true));
		}

	}

	~XMatrixTrees() = default;

	void Update(const Wavefunction& Psi, const Tree& tree) {
		using namespace TreeFunctions;
		assert(xops_.size() == xmats_.size());
		assert(xops_.size() == xholes_.size());
		Represent(xmats_, xops_, Psi, Psi, tree);
		Contraction(xholes_, xmats_, Psi, Psi, tree);
	}

	void print() const {
		cout << "Xs:\n";
		for (const auto& x : xmats_) {
			x.print();
		}

		for (const auto& xhole : xholes_) {
			xhole.print();
		}
	}

	size_t size() const {
		size_t n = xops_.size();
		assert(n == xmats_.size());
		assert(n == xholes_.size());
		return n;
	}

	SOPcd xops_;
	SparseMatrixTreescd xmats_;
	SparseMatrixTreescd xholes_;
};

#endif //XMATRIXTREES_H
