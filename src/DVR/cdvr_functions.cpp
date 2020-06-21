//
// Created by Roman Ellerbrock on 3/18/20.
//

#include "cdvr_functions.h"

namespace cdvr_functions {

	/**
	 * Notes on notation:
	 * Cdown: Top-down normalized, reads C^{p(p\circ k)} in equations
	 * Cup: Bottom-up normalized, reads C^{p(p\circ 0)} in equations
	 */

	/**
	 * Fill a vector with grid points corresponding to intex "idx".
	 * @param X
	 * @param idx
	 * @param grids
	 * @param node
	 */

	void fillX(Vectord& X, size_t idx, const TreeGrids& grids, const Node& node) {
		/// Fill full-dimensional grid points from all neighboring node-grids
		for (size_t k = 0; k < X.Dim(); ++k) {
			const SparseVectorTreed& grid = grids[k];
			if (grid.Active(node)) {
				const Vectord& localgrid = grid[node];
				X(k) = localgrid(idx);
			}
		}
	}

	void fillXNode(Vectord& X, vector<size_t> idx, const TreeGrids& grids,
		const TreeGrids& holegrids, const Node& node) {
		/// \brief Fill X for a node-grid
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

	void fillXEdge(Vectord& X, vector<size_t> idx, const TreeGrids& grids,
		const TreeGrids& holegrids, const Node& node) {
		/// \brief Fill X for an edge-grid
		assert(grids.size() == holegrids.size());
		assert(X.Dim() == grids.size());
		assert(idx.size() == 2);
		fillX(X, idx.front(), grids, node);
		fillX(X, idx.back(), holegrids, node);
	}

	void CalculateDeltaEdgeBottom(Tensorcd& deltaV, const Tensorcd& Cup,
		const Tensorcd& Vnode, const Matrixd& Vedge, const Node& node) {
		/**
		 * \brief Calculate deltaV-tensor for bottomlayer node
		 * @param deltaV "matrix"-representation of operator
		 * @param C Wavefunction coefficients at node
		 * @param Vnode Potential evaluated at node-DVR
		 * @param Vedge Potential evaluated at edge-DVR
		 * @param node Node at which the edge is pointing to (std convention for edges in trees)
		 *
		 * There are no underlying deltaV matrices at node. Therefore, a simplified expression.
		 */

		const TensorShape& shape = Vnode.shape();
		const TensorShape& Cshape = Cup.shape();
		size_t dim = shape.lastDimension();
		size_t dimbef = shape.lastBefore();
		assert(shape.totalDimension() == Cshape.totalDimension());

		/// Sanity checks
		assert(!node.isToplayer());
		assert(node.isBottomlayer());
		assert(deltaV.shape().totalDimension() == pow(shape.lastDimension(), 4));

		deltaV.Zero();
		for (size_t l1 = 0; l1 < dim; ++l1) {
			for (size_t i0 = 0; i0 < dim; ++i0) {
				for (size_t m1 = 0; m1 < dim; ++m1) {
					vector<size_t> idxs = {l1, i0, m1, i0};
					size_t J = indexMapping(idxs, deltaV.shape());
					for (size_t Ibef = 0; Ibef < dimbef; Ibef++) {
						deltaV(J) += conj(Cup(Ibef, l1)) * Vnode(Ibef, i0) * Cup(Ibef, m1);
					}
				}
			}
		}

		/// Substract Vedge from diagonal
		for (size_t l0 = 0; l0 < dim; ++l0) {
			for (size_t l1 = 0; l1 < dim; ++l1) {
				vector<size_t> idxs = {l0, l1, l0, l1};
				size_t L = indexMapping(idxs, deltaV.shape());
				deltaV(L) -= Vedge(l0, l1);
			}
		}
	}

	void CalculateDeltaVs(DeltaVTree& deltaVs, const ExplicitEdgeWavefunction& Chi,
		const TensorTreecd& Vnodes, const MatrixTreed& Vedges,
		const Tree& tree) {

		const TensorTreecd& Cup = Chi.BottomUpNormalized(tree);
		for (const Node& node : tree) {
			if (node.isBottomlayer()) {
				CalculateDeltaEdgeBottom(deltaVs[node], Cup[node],
					Vnodes[node], Vedges[node], node);
			} else if (!node.isToplayer()) {
			}
		}
	}

	Tensorcd ApplyCorrection(const Tensorcd& Phi, const Tensorcd& C,
		const Tensorcd& deltaV, const Node& child) {

		size_t k = child.childIdx();
		const TensorShape& shape = deltaV.shape();

		child.info();
		/// I: Contract over
		Matrixcd X = Contraction(C, Phi, k);
		cout << "X:\n";
		X.print();

		/// II: Apply DeltaV
		size_t dim = shape[k];
		Matrixcd Y(dim, dim);
		for (size_t L = 0; L < shape.totalDimension(); ++L) {
			const auto l = indexMapping(L, shape);
			Y(l[0], l[1]) += deltaV(L) * X(l[2], l[3]);
		}
		cout << "Y:\n";
		Y.print();

		/// III: M * C
		Tensorcd VPhi = MatrixTensor(Y, C, k);
		cout << "VPhi:\n";
		VPhi.print();
		return VPhi;
	}

	Tensorcd Apply(const Tensorcd& Phi, const Tensorcd& V,
		const TensorTreecd& Cdown, const DeltaVTree& deltaVs, const Node& node) {

		Tensorcd VPhi = coeffprod(V, Phi);

		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				VPhi += ApplyCorrection(Phi, Cdown[child], deltaVs[child], child);
			}
		}
		return VPhi;
	}

	TensorTreecd Apply(const ExplicitEdgeWavefunction& Chi,
		const TensorTreecd& V, const DeltaVTree& DeltaVs, const Tree& tree) {

		const TensorTreecd& C = Chi.TopDownNormalized(tree);
		TensorTreecd Psi = Chi.BottomUpNormalized(tree);
		TensorTreecd VPsi = Psi;

		for (const Node& node : tree) {
			VPsi[node] = Apply(Psi[node], V[node], C, DeltaVs, node);
		}
		return VPsi;
	}
}

