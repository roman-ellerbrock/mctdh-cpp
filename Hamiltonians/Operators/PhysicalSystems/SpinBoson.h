//
// Created by Roman Ellerbrock on 8/9/20.
//

#ifndef SPINBOSON_H
#define SPINBOSON_H
#include "TreeOperators/LeafFunction.h"
#include "TreeOperators/SumOfProductsOperator.h"

namespace Operator {

	SOPcd Exciton(const string& filename, const Tree& tree);

	class PrimitiveProjector: public LeafOperatorcd {
	public:
		explicit PrimitiveProjector(size_t idx)	: idx_(idx) {}
		~PrimitiveProjector() = default;
		void Apply(const LeafInterface& grid, Tensorcd& PPhi, const Tensorcd& Phi) const override {
			const TensorShape& tdim = Phi.shape();
			PPhi.Zero();
			assert(idx_ < tdim.lastBefore());
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				PPhi(idx_, n) = Phi(idx_, n);
			}
		}
	private:
		size_t idx_;
	};

	class VectorProjector: public LeafOperatorcd {
	public:
		explicit VectorProjector(const Vectord& U) : U_(U) {
		}

		void Apply(const LeafInterface& grid, Tensorcd& PPhi, const Tensorcd& Phi) const override {
			const TensorShape& tdim = Phi.shape();
			assert(idx_ < tdim.lastBefore());
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				for (size_t i = 0; i < tdim.lastBefore(); ++i) {
					PPhi(i, n) = U_(i) * Phi(i, n);
				}
			}
		}

	protected:
		Vectord U_;
	};

}

#endif //SPINBOSON_H
