//
// Created by Roman Ellerbrock on 4/3/20.
//

#include "ExplicitEdgeWavefunction.h"
#include "TreeClasses/SpectralDecompositionTree.h"

ExplicitEdgeWavefunction::ExplicitEdgeWavefunction(const Wavefunction& Psi, const Tree& tree) {
	Wavefunction& nodes_ = first;
	MatrixTreecd& edges_ = second;
	nodes_ = Psi;

	/// Transform to canonical representation
	CanonicalTransformation(nodes_, tree, true);

	/// Get edge matrices
	MatrixTreecd rho_tmp = TreeFunctions::Contraction(nodes_, tree, true);
	SpectralDecompositionTreecd rho_x(rho_tmp, tree);
	auto sqrt_rho_x = Sqrt(rho_x);

	edges_ = to_matrixtree(sqrt_rho_x, tree);

	/// counter-rotate nodes
	MatrixTreecd inv_sqrt_rho = sqrt_rho_x.Invert(tree);
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const Node& parent = node.parent();
			size_t childidx = node.childIdx();
			nodes_[parent] = MatrixTensor(inv_sqrt_rho[node], nodes_[parent], childidx);
		}
	}

	/// Check contractions
	cout << "nodes:\n";
	for (const Node& node : tree) {
		if (!node.isBottomlayer()) {
			node.info();
			const auto& phi = nodes_[node];
			size_t d = node.nChildren();
			for (size_t i = 0; i < d; ++i) {
				auto xphi(phi);
				for (size_t k = 0; k < d; ++k) {
					const Node& child = node.child(k);
					if (k != i) {
						xphi = MatrixTensor(edges_[child], xphi, k);
					}
				}
				if (!node.isToplayer()) {
					const Node& parent = node.parent();
					xphi = MatrixTensor(edges_[parent], xphi, d);
				}
				cout << "contraction " << i << endl;
				auto rho = Contraction(xphi, xphi, i);
				rho.print();
			}
		}
	}

}
