//
// Created by Roman Ellerbrock on 3/8/20.
//
#include "Core/tLayerInterface.h"

template <typename T>
double tLayerInterface<T>::Error(const Tensor<T>& Phi, const Tensor<T>& Chi)const {
	double Delta = 0;
	if (!node_->isToplayer()) {
		// Density weighted Error
		const auto& rho = hRep_->rho_;
		const auto& rhomat = rho[*node_];
		double norm = real(rhomat.trace());
/*		auto C = Phi - Chi;
		C = multStateAB(rhomat, C);
		for (size_t j = 0; j < C.shape().totalDimension(); j++)
			Delta += pow(abs(C(j)), 2);
		Delta /= norm;
		*/
		size_t k = Phi.shape().lastIdx();
		auto pPhi = matrixTensor(rhomat, Phi, k);
		auto pChi = matrixTensor(rhomat, Phi, k);
		auto scc = Chi.dotProduct(pChi);
		auto spp = Phi.dotProduct(pPhi);
		auto spc = Phi.dotProduct(pChi);
		auto scp = Chi.dotProduct(pPhi);
		Matrix<T> s = scc + spp - spc - scp;
		Delta = s.frobeniusNorm() / norm;
	} else {
		// @TODO: Check if this is ok:
		// Incorporating norm
		auto S = Phi.dotProduct(Phi);
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

template class tLayerInterface<double>;
