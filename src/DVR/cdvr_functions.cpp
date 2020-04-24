//
// Created by Roman Ellerbrock on 3/18/20.
//

#include "cdvr_functions.h"

namespace cdvr_functions {

	void CalculateLocal(Tensorcd& deltaV, const Tensorcd& Cup,
		const Tensorcd& Cdown, const Tensorcd& V_dvr,
		const Matrixd& V_edge, const Node& node) {

		const TensorShape& shape = V_dvr.shape();

		deltaV.Zero();
		assert(deltaV.shape().totalDimension() == pow(shape.lastDimension(), 3));

		for (size_t j = 0; j < shape.lastDimension(); ++j) {
			for (size_t l = 0; l < shape.lastDimension(); ++l) {
				for (size_t m = 0; m < shape.lastDimension(); ++m) {

					/// Add simple grid
					vector<size_t> idxs = {m, l, j};
					size_t J = indexMapping(idxs, deltaV.shape());

					for (size_t I = 0; I < shape.lastBefore(); ++I) {
						deltaV(J) += conj(Phi(I, m)) * V_dvr(I, j) * Phi(I, l);
					}

				}

				/// substract edge grid
				vector<size_t> idxs = {l, l, j};
				size_t J = indexMapping(idxs, deltaV.shape());
				deltaV(J) -= V_edge(l, j);
			}
		}

		/// @TODO: Apply underlying deltaVs

	}

	void Update(DeltaVTree& deltaVs, const ExplicitEdgeWavefunction& Chi,
		const TensorTreecd& V_dvr, const MatrixTreed& V_edge_dvr,
		const Tree& tree) {

		for (const Node& node : tree) {
			if (!node.isToplayer()) {
				CalculateLocal(deltaVs[node], Psi[node], V_dvr[node], V_edge_dvr[node], node);
			}
		}
	}

	void ApplyCorrection(Tensorcd& VPhi, const Tensorcd& A, const Tensorcd& deltaV,
		const Node& child) {

		const TensorShape& shape = A.shape();
		size_t k = child.childIdx();

		for (size_t j = 0; j < shape.lastDimension(); ++j) {

			for (size_t Mbef = 0; Mbef < shape.before(k); ++Mbef) {
				for (size_t Maft = 0; Maft < shape.after(k); ++Maft) {
					for (size_t m = 0; m < shape[k]; ++m) {

						for (size_t Nbef = 0; Nbef < shape.before(k); ++Nbef) {
							for (size_t Naft = 0; Naft < shape.after(k); ++Naft) {
								for (size_t n = 0; n < shape[k]; ++k) {

									size_t delta_idx = indexMapping({m, n, j}, deltaV.shape());
									VPhi(Mbef, m, Maft, k) += A(Mbef, m, Maft, k) *
										deltaV(delta_idx) * conj(A(Nbef, n, Naft, k)) * A(Nbef, n, Naft, k);
								}
							}
						}
					}
				}
			}
		}
	}

	void Apply(Tensorcd& VPhi, const Tensorcd& Phi, const Tensorcd& V,
		const DeltaVTree& deltaVs, const Node& node) {

		const TensorShape& shape = node.shape();
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			VPhi(I) += V(I) * Phi(I);
		}

		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				ApplyCorrection(VPhi, Phi, deltaVs[child], child);
			}
		}
	}
}