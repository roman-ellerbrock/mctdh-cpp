//
// Created by Roman Ellerbrock on 10/31/22.
//

#ifndef QUADRATICPARTITIONS_H
#define QUADRATICPARTITIONS_H
#include "TreeClasses/MatrixTreeFunctions.h"

class QuadraticPartitions : public NodeAttribute<vector<size_t>>{
public:

	QuadraticPartitions() = default;
	~QuadraticPartitions() = default;

	void initialize(const Tree& tree) {
		attributes_.clear();
		for (const Node& node : tree) {
			attributes_.push_back(vector<size_t>());
		}
	}

	QuadraticPartitions(const Tree& tree) {
		initialize(tree);
	}

	void build(const Tree& tree) {
		for (const Node& node : tree) {
			vector<size_t>& nodeidx = operator[](node);
			if (node.isBottomlayer()) {
				const Leaf& leaf = node.getLeaf();
				size_t mode = leaf.mode();
				nodeidx.push_back(mode);
			} else {
				/// add indices of every child
				for (size_t k = 0; k < node.nChildren(); ++k) {
					const Node& child = node.child(k);
					const vector<size_t> childidx = operator[](child);
					nodeidx.insert(nodeidx.end(), childidx.begin(), childidx.end());
				}
			}
		}
	}

	void print(const Tree& tree) const {
		for (const Node& node : tree) {
			node.info();
			for (size_t x : operator[](node)) {
				cout << x << "\t";
			}
			cout << endl;
		}
	}

};

#endif //QUADRATICPARTITIONS_H
