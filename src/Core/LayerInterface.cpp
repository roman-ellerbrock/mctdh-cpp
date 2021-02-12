//
// Created by Roman Ellerbrock on 3/8/20.
//
#include "Core/LayerInterface.h"

double LayerInterface::Error(const Tensorcd& Phi, const Tensorcd& Chi)const {
	double Delta = 0;
	if (!node_->isToplayer()) {
		// Density weighted Error
		const MatrixTreecd& rho = hRep_->rho_;
		const Matrixcd& rhomat = rho[*node_];
		double norm = real(rhomat.trace());
		Tensorcd C = Phi - Chi;
		C = multStateAB(rhomat, C);
		for (size_t j = 0; j < C.shape().totalDimension(); j++)
			Delta += pow(abs(C(j)), 2);
		Delta /= norm;
	} else {
		// @TODO: Check if this is ok:
		// Incorporating norm
		Matrixcd S = Phi.dotProduct(Phi);
		double norm = real(S.trace());
		for (size_t j = 0; j < Phi.shape().totalDimension(); j++) {
			// Primitive Error without density matrix stuff
			Delta += pow(abs(Phi(j) - Chi(j)), 2);
		}
		norm /= Phi.shape().lastDimension();
		Delta /= norm;
	}
	return sqrt(Delta);
}

