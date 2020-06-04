//
// Created by Roman Ellerbrock on 3/13/20.
//

#include "CDVR.h"

CDVR::CDVR(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part)
	: tddvr_(Psi, tree), nodedvr_(tree), edgedvr_(tree), deltaV_(tree) {
	Update(Psi, V, tree, part);
}

void fillX(Vectord& X, size_t idx, const TreeGrids& grids, const Node& node) {
	for (size_t k = 0; k < X.Dim(); ++k) {
		const SparseVectorTreed& grid = grids[k];
		if (grid.Active(node)) {
			const Vectord& localgrid = grid[node];
			X(k) = localgrid(idx);
		}
	}
}

void fillXdvr(Vectord& X, vector<size_t> idx, const TreeGrids& grids, const TreeGrids& holegrids,
	const Node& node) {
	assert(grids.size() == holegrids.size());
	assert(X.Dim() == grids.size());
	if (node.isBottomlayer()) {
		const Leaf& leaf = node.getLeaf();
		const auto& g = leaf.PrimitiveGrid();
		assert(g.HasDVR());
		const Vectord& x = g.GetX();
		X(leaf.Mode()) = x(idx.front());
	} else {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			fillX(X, idx[k], grids, child);
		}
	}

	fillX(X, idx.back(), holegrids, node);
}

void fillXcorrection(Vectord& X, vector<size_t> idx, const TreeGrids& grids, const TreeGrids& holegrids,
	const Node& node) {
	assert(grids.size() == holegrids.size());
	assert(X.Dim() == grids.size());
	assert(idx.size() == 2);
	fillX(X, idx.front(), grids, node);
	fillX(X, idx.back(), holegrids, node);
}

void UpdateDVRLocal(Tensorcd& dvr, const TreeGrids& grids, const TreeGrids& holegrids,
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
		fillXdvr(X, idxs, grids, holegrids, node);
		dvr(I) = V.Evaluate(X, part);
	}
}

void UpdateDVR(TensorTreecd& dvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Tree& tree, size_t part) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			UpdateDVRLocal(dvr[node], grids, holegrids, V, node, part);
		}
	}
}

void UpdateCDVRLocal(Matrixd& cdvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Node& node, size_t part) {
	const TensorShape& shape = node.shape();
	size_t ngrid = shape.lastDimension();
	Vectord X(grids.size());
	TensorShape gridshape({ngrid, ngrid});
	for (size_t I = 0; I < gridshape.totalDimension(); ++I) {
		auto idxs = indexMapping(I, gridshape);
		fillXcorrection(X, idxs, grids, holegrids, node);
		cdvr(idxs.front(), idxs.back()) = V.Evaluate(X, part);
	}
}

void UpdateCDVR(MatrixTreed& cdvr, const TreeGrids& grids, const TreeGrids& holegrids,
	const Potential& V, const Tree& tree, size_t part) {

	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			UpdateCDVRLocal(cdvr[node], grids, holegrids, V, node, part);
		}
	}
}

void CDVR::Update(const Wavefunction& Psi, const Potential& V, const Tree& tree, size_t part) {

	/// Build X-matrices, diagonalize them simultaneously
	tddvr_.Update(Psi, tree);

	UpdateDVR(nodedvr_, tddvr_.grids_, tddvr_.hole_grids_, V, tree, part);

	UpdateCDVR(edgedvr_, tddvr_.grids_, tddvr_.hole_grids_, V, tree, part);

//	cdvr_functions::Update(deltaV_, Psi, dvr_, cdvr_, tree);

}

