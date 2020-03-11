//
// Created by Roman Ellerbrock on 3/11/20.
//
#include "DVR/XMatrixTrees.h"

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
