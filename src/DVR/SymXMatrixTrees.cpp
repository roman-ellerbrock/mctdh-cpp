//
//

#include "SymXMatrixTrees.h"

Tensorcd SymXMatrixTrees::OptimizeUp(const Tensorcd& Phi, const Matrixcd& rho,
	const Node& node, const Node& node_small) const {

	size_t n_occ = node_small.shape().lastDimension();

/*	auto X = BuildX(Phi, rho, node);
	X = UnProject(n_occ, X, Phi);
	auto xspec = Diagonalize(X);
	auto oPhi = Occupy(Phi, xspec.first, n_occ, node);

	return oPhi;
 */
}

MatrixTensorTree SymXMatrixTrees::OptimizeUp(const MatrixTensorTree& Psi,
	const Tree& tree, const Tree& tree_small) const {
	TensorTreecd PsiUp = Psi.BottomUpNormalized(tree);
	MatrixTreecd rho = TreeFunctions::contraction(PsiUp, tree, true);
	for (const Node& node : tree) {
		const Node& node_small = tree_small.getNode(node.address());
		PsiUp[node] = OptimizeUp(PsiUp[node], rho[node], node, node_small);
	}
	return MatrixTensorTree(PsiUp, tree, true);
}
