//
// Created by Roman on 3/31/2019.
//

#include "CoupledHO.h"

SOPcd CoupledHO(const Tree& tree) {
	LeafFuncd x = &LeafInterface::applyX;
	LeafFuncd x2 = &LeafInterface::applyX2;
	LeafFuncd kin = &LeafInterface::applyKin;
	LeafFuncd p = &LeafInterface::applyP;

	constexpr double cm = 219474.6313705;
	constexpr double lambda = 2000. / cm;
	constexpr double omega = 4000. / cm;

	constexpr double c = 0.5 * omega * omega;

	size_t f = tree.nLeaves();

	SOPcd H;
// Kinetic energy
	for (size_t k = 0; k < f; ++k) {
		{
			MLOcd M;
			M.push_back(kin, k);
			H.push_back(M, 1.);
		}
		{
			MLOcd M;
			M.push_back(x2, k);
			H.push_back(M, c);
		}
	}

	for (size_t k = 0; k < f; ++k) {
		size_t kn = (k + 1) % f;
		MLOcd M;
		M.push_back(x, k);
		M.push_back(x, kn);
		H.push_back(M, lambda * lambda);
	}

	return H;

	/// H = a * (h1 * h2 *...*hd) + b(g1 * g2 * ...)
}



