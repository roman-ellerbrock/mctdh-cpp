//
// Created by Roman Ellerbrock on 3/18/20.
//

#include "CDVREdgeTensor.h"

void CDVREdgeTensor::Initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const TensorShape& shape = node.shape();
			size_t dim = shape.lastDimension();
			TensorShape new_shape({dim, dim, dim});
			attributes_.emplace_back(new_shape);
		}
	}
}

void CalculateLocal(Tensorcd& edge_dvr, const Tensorcd& Phi,
	const Tensorcd& V_dvr, const MatrixTreed& V_edge_dvr,
	const Node& node) {

	const Matrixd& V_edge = V_edge_dvr[node];
	const TensorShape& shape = V_dvr.shape();

	edge_dvr.Zero();

	for (size_t j = 0; j < shape.lastDimension(); ++j) {
		for (size_t m = 0; m < shape.lastDimension(); ++m) {
			for (size_t l = 0; l < shape.lastDimension(); ++l) {

				vector<size_t> idxs = {j, l, m};
				size_t J = indexMapping(idxs, shape);

				for (size_t I = 0; I < shape.lastBefore(); ++I) {
					edge_dvr(J) += conj(Phi(I, m)) * V_dvr(I, j) * Phi(I, l);
				}

			}
			vector<size_t> idxs = {j, m, m};
			size_t J = indexMapping(idxs, shape);
			edge_dvr(J) -= V_edge(m, j);
		}
	}

}

void CDVREdgeTensor::Calculate(const Wavefunction& Psi,
	const TensorTreecd& V_dvr, const MatrixTreed& V_edge_dvr,
	const Tree& tree) {

	for (const Node& node : tree) {
		if(!node.isToplayer()) {
			CalculateLocal((*this)[node], Psi[node], V_dvr[node], V_edge_dvr, node);
		}
	}

}
