//
// Created by Roman Ellerbrock on 10/31/22.
//

#ifndef QUADRATICHOLEPARTITIONS_H
#define QUADRATICHOLEPARTITIONS_H
#include "QuadraticPartitions.h"

class QuadraticHolePartitions : public NodeAttribute<vector<size_t>> {
public:
	QuadraticHolePartitions() = default;
	~QuadraticHolePartitions() = default;

	QuadraticHolePartitions(const Tree& tree) {
		initialize(tree);
	}

	void initialize(const Tree& tree) {
		attributes_.clear();
		attributes_.resize(tree.nNodes());
	}

	void build(const QuadraticPartitions& idx, const Tree& tree) {
		/// skip root and go top-down
		for (int i = tree.nNodes() - 2; i >= 0; i--) {
			const Node& node = tree.nodes()[i];
			const Node& parent = node.parent();
			vector<size_t>& nodeidx = operator[](node);
			for (size_t k = 0; k < parent.nChildren(); ++k) {
				if (k == node.childIdx()) { continue; }
				const Node& child = parent.child(k);
				const vector<size_t>& parentchildidx = idx[child];
				nodeidx.insert(nodeidx.end(), parentchildidx.begin(), parentchildidx.end());
			}
			if (!parent.isToplayer()) {
				const vector<size_t>& parentidx = operator[](parent);
				nodeidx.insert(nodeidx.end(), parentidx.begin(), parentidx.end());
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


#endif //QUADRATICHOLEPARTITIONS_H
