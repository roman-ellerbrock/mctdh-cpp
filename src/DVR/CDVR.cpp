//
// Created by Roman Ellerbrock on 3/13/20.
//

#include "CDVR.h"

using namespace cdvr_functions;

/*CDVR::CDVR(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree, size_t part)
	: ltree_(tree), tddvr_(Psi, tree), Vnode_(tree), Vedge_(tree),
		deltaV_(tree), Chi_(Psi, tree, true) {
	auto Chi = Psi;
	TreeFunctions::Adjust(Chi, ltree_);
	Update(Chi, V, ltree_, part);
}*/

CDVR::CDVR(const Tree& tree)
	: ltree_(tree), tddvr_(tree), Vnode_(tree),
	  Vedge_(tree), deltaV_(tree) {
}

void UpdateNodeDVRLocal(Tensorcd& dvr, const TreeGrids& grids,
	const TreeGrids& holegrids, const PotentialOperator& V,
	const Node& node, size_t part, bool out = false, ostream& os = cout) {

	/// Check-a-lot
	assert(grids.size() == holegrids.size());
	const TensorShape& shape = dvr.shape();
	for (size_t k = 0; k < grids.size(); ++k) {
		const SparseVectorTreed& grid = grids[k];
		const SparseVectorTreed& holegrid = holegrids[k];
		assert(grid.Active(node) != holegrid.Active(node));
	}

	Vectord X(grids.size());
	vector<size_t> idxs(shape.order());
	for (size_t I = 0; I < shape.totalDimension(); ++I) {
//		auto idxs = indexMapping(I, shape);
		indexMapping(idxs, I, shape);
		fillXNode(X, idxs, grids, holegrids, node);
		dvr(I) = V.Evaluate(X, part);

		if (out) {
			for (size_t i = 0; i < X.Dim(); ++i) {
				os << X(i) << "\t";
			}
			os << real(dvr(I)) << endl;
		}
	}
}

void UpdateNodeDVR(TensorTreecd& dvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const PotentialOperator& V, const Tree& tree, size_t part, bool out, ostream& os) {
	for (const Node& node : tree) {
		UpdateNodeDVRLocal(dvr[node], grids, holegrids, V, node, part, out, os);
	}
}

void UpdateEdgeDVRLocal(Matrixd& edgedvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const PotentialOperator& V, const Node& node, size_t part, bool out = false, ostream& os = cout) {
	const TensorShape& shape = node.shape();
	size_t ngrid = shape.lastDimension();
	Vectord X(grids.size());
	TensorShape gridshape({ngrid, ngrid});
	vector<size_t> idxs(gridshape.order());
	for (size_t I = 0; I < gridshape.totalDimension(); ++I) {
//		auto idxs = indexMapping(I, gridshape);
		indexMapping(idxs, I, gridshape);
		fillXEdge(X, idxs, grids, holegrids, node);
		edgedvr(idxs.front(), idxs.back()) = V.Evaluate(X, part);
		if (out) {
			for (size_t i = 0; i < X.Dim(); ++i) {
				os << X(i) << "\t";
			}
			os << edgedvr(idxs.front(), idxs.back()) << endl;
		}
	}
}

void UpdateEdgeDVR(MatrixTreed& cdvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const PotentialOperator& V, const Tree& tree, size_t part, bool out, ostream& os = cout) {

	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			UpdateEdgeDVRLocal(cdvr[node], grids, holegrids, V, node, part, out, os);
		}
	}
}

void CDVR::Update(const Wavefunction& Psi, const PotentialOperator& V,
	const Tree& tree, size_t part, bool out, ostream& os) {

	/// Get Edge wavefunction
	Chi_ = ExplicitEdgeWavefunction(Psi, tree, true);

	/// Build X-matrices, diagonalize them simultaneously
	tddvr_.Update(Psi, tree);

	/// Transform to grid
	tddvr_.GridTransformation(Chi_, tree);

	/// Save top-down normalized wavefunction, since it is needed to apply the CDVR-operator
	Cdown_ = Chi_.TopDownNormalized(tree);

	/// Evaluate potential at Nodes and edges
	UpdateNodeDVR(Vnode_, tddvr_.grids_, tddvr_.hole_grids_, V, tree, part, out, os);
	UpdateEdgeDVR(Vedge_, tddvr_.grids_, tddvr_.hole_grids_, V, tree, part, out, os);
	if (out) { os << endl; }

	/// Evaluate correction matrices
	cdvr_functions::CalculateDeltaVs(deltaV_, Chi_, Vnode_, Vedge_, tree);
}

//Tensorcd CDVR::Apply(Tensorcd Phi, const Matrixcd& sqrho, const Node& node) const {
Tensorcd CDVR::Apply(Tensorcd Phi, const SpectralDecompositioncd& rho_x, const Node& node) const {

	auto x_sqrt_rho = sqrt(rho_x);
	x_sqrt_rho.second = Regularize(x_sqrt_rho.second, 1e-6);
	auto sqrho = toMatrix(x_sqrt_rho);

	if (!node.isToplayer()) {
		Phi = MatrixTensor(sqrho, Phi, node.nChildren());
	}
	/// adjust size
	const Node& lnode = ltree_.GetNode(node.Address());
	if (lnode.shape() != Phi.shape()) {
		Phi = Phi.AdjustDimensions(lnode.shape());
	}

	tddvr_.NodeTransformation(Phi, node, false);

	auto VXi = cdvr_functions::Apply(Phi, Vnode_[node], Cdown_, deltaV_, node);

	tddvr_.NodeTransformation(VXi, node, true);

	/// adjust size
	if (node.shape() != VXi.shape()) {
		VXi = VXi.AdjustDimensions(node.shape());
	}

	if (!node.isToplayer()) {
//		VXi = TensorMatrix(VXi, sqrho, node.nChildren());
		VXi = MatrixTensor(sqrho.Adjoint(), VXi, node.nChildren());
	}
	return VXi;
}

void CDVR::Update2(Wavefunction Psi, const PotentialOperator& V,
	const Tree& smalltree, size_t part, bool out, ostream& os) {

	/// Inflate wavefunction basis
	TreeFunctions::Adjust(Psi, ltree_);
//	tddvr_.Update(Psi, ltree_);
	tddvr_.Xs_.Update(Psi, ltree_);
	auto rho = TreeFunctions::Contraction(Psi, ltree_, true);
	Psi = tddvr_.Xs_.Optimize(Psi, rho, ltree_, smalltree);

	/// Build X-matrices, diagonalize them simultaneously
	tddvr_.Update(Psi, ltree_);

	/// Get Edge wavefunction
	Chi_ = ExplicitEdgeWavefunction(Psi, ltree_, true);

	/// Transform to grid
	tddvr_.GridTransformation(Chi_, ltree_);

	/// Save top-down normalized wavefunction, since it is needed to apply the CDVR-operator
	Cdown_ = Chi_.TopDownNormalized(ltree_);

	/// Evaluate potential at Nodes and edges
	UpdateNodeDVR(Vnode_, tddvr_.grids_, tddvr_.hole_grids_, V, ltree_, part, out, os);
	UpdateEdgeDVR(Vedge_, tddvr_.grids_, tddvr_.hole_grids_, V, ltree_, part, out, os);
	if (out) { os << endl; }

	/// Evaluate correction matrices
	cdvr_functions::CalculateDeltaVs(deltaV_, Chi_, Vnode_, Vedge_, ltree_);
}


