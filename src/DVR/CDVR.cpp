//
// Created by Roman Ellerbrock on 3/13/20.
//

#include "CDVR.h"

using namespace cdvr_functions;

CDVR::CDVR(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part)
	: tddvr_(Psi, tree), Vnode_(tree), Vedge_(tree), deltaV_(tree), Chi_(Psi, tree, true) {
	Update(Psi, V, tree, part);
}

void UpdateNodeDVRLocal(Tensorcd& dvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Node& node, size_t part) {

	/// Check-a-lot
	assert(grids.size() == holegrids.size());
	const TensorShape& shape = dvr.shape();
	for (size_t k = 0; k < grids.size(); ++k) {
		const SparseVectorTreed& grid = grids[k];
		const SparseVectorTreed& holegrid = holegrids[k];
		assert(grid.Active(node) != holegrid.Active(node));
	}

	Vectord X(grids.size());
	for (size_t I = 0; I < shape.totalDimension(); ++I) {
		auto idxs = indexMapping(I, shape);
		fillXNode(X, idxs, grids, holegrids, node);
		dvr(I) = V.Evaluate(X, part);
	}
}

void UpdateNodeDVR(TensorTreecd& dvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Tree& tree, size_t part) {
	for (const Node& node : tree) {
		UpdateNodeDVRLocal(dvr[node], grids, holegrids, V, node, part);
	}
}

void UpdateEdgeDVRLocal(Matrixd& edgedvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Node& node, size_t part) {
	const TensorShape& shape = node.shape();
	size_t ngrid = shape.lastDimension();
	Vectord X(grids.size());
	TensorShape gridshape({ngrid, ngrid});
	for (size_t I = 0; I < gridshape.totalDimension(); ++I) {
		auto idxs = indexMapping(I, gridshape);
		fillXEdge(X, idxs, grids, holegrids, node);
		edgedvr(idxs.front(), idxs.back()) = V.Evaluate(X, part);
	}
}

void UpdateEdgeDVR(MatrixTreed& cdvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Tree& tree, size_t part) {

	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			UpdateEdgeDVRLocal(cdvr[node], grids, holegrids, V, node, part);
		}
	}
}

void CDVR::Update(const Wavefunction& Psi, const Potential& V,
	const Tree& tree, size_t part) {

	/// Get Edge wavefunction
	Chi_ = ExplicitEdgeWavefunction(Psi, tree, true);

	/// Build X-matrices, diagonalize them simultaneously
	tddvr_.Update(Psi, tree);

	/// Transform to grid
	tddvr_.GridTransformation(Chi_, tree);

	/// Save top-down normalized wavefunction, since it is needed to apply the CDVR-operator
	Cdown_ = Chi_.TopDownNormalized(tree);

	/// Evaluate potential at Nodes and edges
	UpdateNodeDVR(Vnode_, tddvr_.grids_, tddvr_.hole_grids_, V, tree, part);
	UpdateEdgeDVR(Vedge_, tddvr_.grids_, tddvr_.hole_grids_, V, tree, part);

	/// Evaluate correction matrices
	cdvr_functions::CalculateDeltaVs(deltaV_, Chi_, Vnode_, Vedge_, tree);
}

Tensorcd CDVR::Apply(Tensorcd Phi, const Node& node) const {

	tddvr_.NodeTransformation(Phi, node, false);

	auto VXi = cdvr_functions::Apply(Phi, Vnode_[node], Cdown_, deltaV_, node);

	tddvr_.NodeTransformation(VXi, node, true);
	return VXi;
}

TensorTreecd CDVR::Apply(const Wavefunction& Psi, const Tree& tree) const {
	ExplicitEdgeWavefunction Chi(Psi, tree, true);

	auto phi = Chi.nodes();
	Wavefunction VPsi2 = Psi;
	for (const Node& node : tree) {
		VPsi2[node] = Apply(phi[node], node);
		tddvr_.NodeTransformation(VPsi2[node], node, false);
	}

	tddvr_.GridTransformation(Chi, tree, false);
	auto VPsi = cdvr_functions::Apply(Chi.nodes(), Cdown_, Vnode_, deltaV_, tree);

	for (const Node& node : tree) {
		const auto xi = Chi.nodes()[node];
		node.info();
		cout << "S1:\n";
		xi.DotProduct(VPsi[node]).print();
		cout << "S2:\n";
		xi.DotProduct(VPsi2[node]).print();
	}

	return VPsi;
}


