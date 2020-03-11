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
		: xops_holes_(Xsop(tree)), xholes_(Xsop(tree), tree) {
		LeafFuncd x = &LeafInterface::applyX;

		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			size_t mode = leaf.Mode();
			MLOcd M(x, mode);

			xmats_.emplace_back(SparseMatrixTreecd(M, tree));
			xops_mats_.push_back(M, 1.);
		}
	}

	~XMatrixTrees() = default;

	void Update(const Wavefunction& Psi, const Tree& tree) {
		using namespace SparseMatrixTreeFunctions;
		assert(xops_mats_.size() == xmats_.size());
		assert(xops_holes_.size() == xholes_.size());
		Represent(xmats_, xops_mats_, Psi, Psi, tree);
		Represent(xholes_, xops_holes_, Psi, Psi, tree);
	}

	void print() const {
		cout << "Xs:\n";
		for (const auto& x : xmats_) {
			x.print();
		}

		xholes_.print();
	}

	size_t size() const {
		size_t n = xops_mats_.size();
		assert(n == xops_holes_.size());
		assert(n == xmats_.size());
		assert(n == xholes_.size());
		return n;
	}

	SOPcd xops_mats_;
	vector<SparseMatrixTreecd> xmats_;

	SOPcd xops_holes_;
	MatrixTreescd xholes_;
};

#endif //XMATRIXTREES_H
