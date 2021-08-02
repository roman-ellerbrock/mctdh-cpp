//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef LAYERINTERFACE_H
#define LAYERINTERFACE_H
#include "Core/tHamiltonianRepresentation.h"

template <typename T>
class tLayerInterface {
public:
	tLayerInterface(const SOP<T>& H, const tHamiltonianRepresentation<T>& hRep,
		const Node& node, T propagation_phase)
		: h_(&H), hRep_(&hRep), node_(&node), propagation_phase_(propagation_phase) {}

	~tLayerInterface() = default;

	void Derivative(double time, Tensor<T>& dPhi, const Tensor<T>& Phi) {
		LayerDerivative(dPhi, time, Phi, *h_, *hRep_, *node_, propagation_phase_);
	}

	double Error(const Tensor<T>& Phi, const Tensor<T>& Chi) const;

	T propagation_phase_;
private:
	const SOP<T>* h_;
	const tHamiltonianRepresentation<T>* hRep_;
	const Node* node_;
};

#endif //LAYERINTERFACE_H
