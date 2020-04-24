#include "mlSOP_model_1.h"

Tensorcd unit(const PrimitiveBasis& grid, const Tensorcd& Acoeffs) {
	return Acoeffs;
}

void mlSOP_model_1::SpecialInitializeBottom(const mctdhBasis& basis)
{
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x2 = &PrimitiveBasis::ApplyX2;
	for (const mctdhNode& node : basis) {
		if (node.IsBottomlayer()) {
			const PhysicalCoordinate& phy = node.PhysCoord();
			size_t mode = phy.Mode();

			MultiParticleOperator M;
			M.push_back(unit, mode);
			M.push_back(x2, mode);
		}
	}
}

void mlSOP_model_1::SpecialInitialize(const mctdhBasis& basis)
{
	for (const mctdhNode& node : basis) {
		const TensorDim& tdim = node.TDim();
		attributes.push_back(Tensorcd(tdim));
	}

	double omega = 2000./219474.6313705;
	double factor = 0.5*omega*omega;

	// factor*x1**2 + factor*x2**2 + ...
	for (const mctdhNode& node : basis) {
		Tensorcd& c = operator[](node);
		if (node.IsBottomlayer()) {
			// 0: 1*unit
			c(0, 0) = 1.;
			// 1: factor*x**2
			c(1, 1) = factor;
		} else if (node.IsToplayer()){ 
			// 0: (x1**2*unit + unit*x2**2) + (x3**2*unit + unit*x4**2)
			c(1, 0) = 1.;
			c(2, 0) = 1.;
		} else {
			// 0: unit*unit
			c(0, 0) = 1.;
			// 1: 0*unit*unit + (x1**2*unit) + (unit*x1**2) + 0*x1**2*x2**2
			c(1, 1) = 1.;
			c(2, 1) = 1.;
		}
	}


}


