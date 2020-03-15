//
// Created by Roman Ellerbrock on 3/9/20.
//

#ifndef TREEGRIDS_H
#define TREEGRIDS_H
#include "TreeClasses/SparseMatrixTreeFunctions.h"

class SparseVectorTreed: public SparseNodeAttribute<Vectord> {
public:
	SparseVectorTreed(const MLOcd& M, const Tree& tree, bool tail = true, bool inverse_tree = false)
		: SparseNodeAttribute<Vectord>(M.targetLeaves(), tree, tail, inverse_tree) {
		SparseVectorTreed::Initialize(tree);
	}

	SparseVectorTreed(const SparseTree& stree, const Tree& tree)
		: SparseNodeAttribute<Vectord>(stree, tree) {
		SparseVectorTreed::Initialize(tree);
	}

	void Initialize(const Tree& tree) override {
		attributes_.clear();
		for (const Node *node_ptr : Active()) {
			const Node& node = *node_ptr;
			const TensorShape& shape = node.shape();
			attributes_.emplace_back(Vectord(shape.lastDimension()));
		}
	}

};

class TreeGrids: public vector<SparseVectorTreed> {
public:
	TreeGrids() = default;

	explicit TreeGrids(const Tree& tree, bool inverse_tree = false) {
		LeafFuncd x = &LeafInterface::applyX;

		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			size_t mode = leaf.Mode();
			MLOcd M(x, mode);
			emplace_back(SparseVectorTreed(M, tree, true, inverse_tree));
		}
	}

	void print(const Tree& tree) const {
		for (const Node& node : tree) {
			node.info();
			for (const auto& grid : *this) {
				if (grid.Active(node)) {
					grid[node].print();
				}
			}
		}
	}
};

#endif //TREEGRIDS_H
