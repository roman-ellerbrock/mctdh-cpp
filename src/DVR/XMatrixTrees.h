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
#include "Core/Tensor_Extension.h"

SOPcd Xsop(const Tree& tree);

Wavefunction Regularize(Wavefunction Psi, const Tree& tree, double eps);

typedef pair<Matrixcd, size_t> MatrixIdx;

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
			holes_.emplace_back(SparseMatrixTreecd(M, tree, false, true));
		}

	}

	~XMatrixTrees() = default;

	void Update(const Wavefunction& Psi, const Tree& tree) {
		using namespace TreeFunctions;
		assert(xops_.size() == mats_.size());
		assert(xops_.size() == holes_.size());
		/**
		 * Analysis:
		 * - density is not included in mean-field x-matrices of lower layers
		 *   in ml trees for inverse trees.
		 * Find Evidence:
		 * - Calculate full mean-field x-matrices (dense tree)
		 *
		 */
		Represent(mats_, xops_, Psi, Psi, tree);
		Wavefunction Chi = Regularize(Psi, tree, 1e-3);
		auto rho = TreeFunctions::Contraction(Chi, tree, true);
		Contraction(holes_, Chi, Chi, mats_, rho, tree);
		UnweightContractions(holes_, Chi, tree);
	}

	void UnweightContractions(vector<SparseMatrixTreecd>& holes,
		const Wavefunction& Psi, const Tree& tree) const {
		auto rho = TreeFunctions::Contraction(Psi, tree, true);
		auto rho_sqrt = sqrt(rho, tree);
		auto isqrt_rho = inverse(rho_sqrt, tree, 1e-5);

		for (auto& xhole : holes) {
			const auto& stree = xhole.Active();
			for (const Node *node_ptr : stree) {
				const Node& node = *node_ptr;
				if (!node.isToplayer()) {
					const auto& isq_rho = isqrt_rho[node];
					auto& x = xhole[node];
//					node.info();
//					x.print();
//					isq_rho.print();
					x = isq_rho * xhole[node] * isq_rho;
//					x.print();
				}
			}
		}
	}

	Tensorcd Optimize(const Tensorcd& Phi, const Matrixcd& rho,
		const Node& node, const Node& node_small) const;

	Wavefunction Optimize(Wavefunction Psi,
		const MatrixTreecd& rho, const Tree& tree, const Tree& tree_small);

	Matrixcd BuildX(const Tensorcd& Phi, const Matrixcd& rho, const Node& node) const;

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

Matrixcd UnProject(size_t n_occupied, const Matrixcd& X,
	const Tensorcd& Phi);

Tensorcd Occupy(const Tensorcd& Phi, const Matrixcd& trafo,
	size_t n_occupied, const Node& node);

#endif //XMATRIXTREES_H
