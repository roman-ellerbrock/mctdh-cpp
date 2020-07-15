//
// Created by Roman Ellerbrock on 3/11/20.
//

#ifndef XMATRIXTREES_H
#define XMATRIXTREES_H
#include "Core/Wavefunction.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/SOPMatrixTrees.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SpectralDecompositionTree.h"

SOPcd Xsop(const Tree& tree);

Wavefunction Regularize(Wavefunction Psi, const Tree& tree, double eps);

class XMatrixTrees {
public:
	explicit XMatrixTrees(const Tree& tree)
		: xops_(Xsop(tree)) {
		LeafFuncd x = &LeafInterface::applyX;

		for (size_t l = 0; l < tree.nLeaves(); ++l) {
			const Leaf& leaf = tree.GetLeaf(l);
			size_t mode = leaf.Mode();
			MLOcd M(x, mode);

			mats_.emplace_back(SparseMatrixTreecd(M, tree));
			holes_.emplace_back(SparseMatrixTreecd(M, tree, true, true));
		}

	}

	~XMatrixTrees() = default;

	void Update(const Wavefunction& Psi, const Tree& tree) {
		using namespace TreeFunctions;
		assert(xops_.size() == mats_.size());
		assert(xops_.size() == holes_.size());
		Represent(mats_, xops_, Psi, Psi, tree);
		Wavefunction Chi = Regularize(Psi, tree, 1e-5);
		Contraction(holes_, mats_, Chi, Chi, tree);
		UnweightContractions(holes_, Chi, tree);
	}

	void UnweightContractions(vector<SparseMatrixTreecd>& holes,
		const Wavefunction& Psi, const Tree& tree) const {
		auto rho = TreeFunctions::Contraction(Psi, tree, true);
		auto rho_sqrt = sqrt(rho, tree);
		auto isqrt_rho = inverse(rho_sqrt, tree, 1e-7);

		for (auto& xhole : holes) {
			const auto& stree = xhole.Active();
			for (const Node *node_ptr : stree) {
				const Node& node = *node_ptr;
				if (!node.isToplayer()) {
					const auto& isq_rho = isqrt_rho[node];
					auto& x = xhole[node];
					x = isq_rho * xhole[node] * isq_rho;
				}
			}
		}
	}

	void print() const {
		cout << "Xs:\n";
		for (const auto& x : mats_) {
			x.print();
		}

		cout << "X holes:\n";
		for (const auto& xhole : holes_) {
			xhole.print();
		}
	}

	size_t size() const {
		size_t n = xops_.size();
		assert(n == mats_.size());
		assert(n == holes_.size());
		return n;
	}

	SOPcd xops_;
	SparseMatrixTreescd mats_;
	SparseMatrixTreescd holes_;
};

#endif //XMATRIXTREES_H
