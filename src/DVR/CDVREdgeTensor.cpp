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

