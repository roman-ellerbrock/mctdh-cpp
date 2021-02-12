//
// Created by Roman Ellerbrock on 3/18/20.
//

#include "cdvr_functions.h"
#include "Core/Tensor_Extension.h"

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
		for (size_t k = 0; k < X.dim(); ++k) {
			const SparseVectorTreed& grid = grids[k];
			if (grid.isActive(node)) {
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
			const auto& g = leaf.interface();
			assert(g.HasDVR());
			const Vectord& x = g.getX();
			X(leaf.mode()) = x(idx.front());
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

	void DeltaEdgeCorrection(Tensorcd& deltaV, const Tensorcd& Cup,
		const TensorTreecd& Cdown, const DeltaVTree& DeltaVs,
		const Node& node) {

		assert(!node.isBottomlayer());
		const TensorShape& shape = node.shape();
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			const Tensorcd& SubDeltaV = DeltaVs[child];
			size_t dimc = shape[k];
			size_t dimn = shape[node.nChildren()];

			auto D = Tensor_Extension::doubleHoleContraction(Cdown[child], Cup, k, node.nChildren());
			auto F = Tensor_Extension::doubleHoleContraction(Cup, Cdown[child], k, node.nChildren());

			TensorShape eshape({dimc, dimc, dimn, dimn});
			Tensorcd E(eshape);
			for (size_t J = 0; J < eshape.totalDimension(); ++J) {
				for (size_t l2 = 0; l2 < dimc; ++l2) {
					for (size_t l1 = 0; l1 < dimc; ++l1) {
						auto idx = indexMapping(J, eshape);
						vector<size_t> vidx({idx[0], idx[1], l1, l2});
						vector<size_t> didx({l2, idx[3], l1, idx[2]});
						E(J) += SubDeltaV(vidx) * D(didx);
					}
				}
			}

			const TensorShape& vshape = deltaV.shape();
			for (size_t I = 0; I < vshape.totalDimension(); ++I) {
				for(size_t l1 = 0; l1 < dimc; ++l1) {
					for (size_t l2 = 0; l2 < dimc; ++l2) {
						auto idx = indexMapping(I, vshape);
						vector<size_t> fidx({l1, idx[0], l2, idx[1]});
						vector<size_t> eidx({l1, l2, idx[2], idx[3]});
						deltaV(I) += F(fidx) * E(eidx);
					}
				}
			}
		}

	}

	void CalculateDeltaEdgeLocal(Tensorcd& deltaV, const Tensorcd& Cup,
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
		assert(deltaV.shape().totalDimension() == pow(shape.lastDimension(), 4));

		deltaV.zero();
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

	void CalculateDeltaVs(DeltaVTree& deltaVs, const MatrixTensorTree& Chi,
		const TensorTreecd& Vnodes, const MatrixTreed& Vedges,
		const Tree& tree) {

		const TensorTreecd& Cup = Chi.BottomUpNormalized(tree);
		const TensorTreecd& Cdown = Chi.TopDownNormalized(tree);
		for (const Node& node : tree) {
			if (node.isBottomlayer()) {
				CalculateDeltaEdgeLocal(deltaVs[node], Cup[node],
					Vnodes[node], Vedges[node], node);

			} else if (!node.isToplayer()) {
				CalculateDeltaEdgeLocal(deltaVs[node], Cup[node],
					Vnodes[node], Vedges[node], node);
				DeltaEdgeCorrection(deltaVs[node], Cup[node], Cdown,
					deltaVs, node);
			}
		}
	}

	Tensorcd ApplyCorrection(const Tensorcd& Phi, const Tensorcd& C,
		const Tensorcd& deltaV, const Node& child) {

		/// I: Contract over
		size_t k = child.childIdx();
		const TensorShape& shape = deltaV.shape();
		assert(C.shape().totalDimension() == Phi.shape().totalDimension());
		size_t dim = Phi.shape()[k];

		Matrixcd X = contraction(C, Phi, k);
		assert(X.Dim1() == dim);
		assert(X.Dim2() == dim);

		/// II: apply DeltaV
		Matrixcd Y(dim, dim);
		for (size_t L = 0; L < shape.totalDimension(); ++L) {
			const auto l = indexMapping(L, shape);
			Y(l[0], l[1]) += deltaV(L) * X(l[3], l[2]);
		}

		/// III: M * C
		return matrixTensor(Y, C, k);
	}

	Tensorcd Apply(const Tensorcd& Xi, const Tensorcd& V,
		const TensorTreecd& Cdown, const DeltaVTree& deltaVs, const Node& node) {

		/// TODO: V diagonal in last idx for toplayer
		Tensorcd VXi = productElementwise(V, Xi);

		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				Tensorcd DeltaXi = ApplyCorrection(Xi, Cdown[child], deltaVs[child], child);
				VXi += DeltaXi;
			}
		}

		return VXi;
	}

	TensorTreecd Apply(const Wavefunction& Psi, const TensorTreecd& Cdown,
		const TensorTreecd& V, const DeltaVTree& DeltaVs, const Tree& tree) {

		Wavefunction VPsi = Psi;
		for (const Node& node : tree) {
			VPsi[node] = Apply(Psi[node], V[node], Cdown, DeltaVs, node);
		}

		return VPsi;
	}
}

