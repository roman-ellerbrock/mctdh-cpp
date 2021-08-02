//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef INTEGRATORINTERFACE_H
#define INTEGRATORINTERFACE_H
#include "Core/tHamiltonianRepresentation.h"

template <typename T>
class IntegratorInterface {
public:
	IntegratorInterface(const SOP<T>& H, const Tree& tree, T phase = 1.)
		: hRep_(H, tree), H_(H), tree_(tree), phase_(phase), dPsi_(tree) {}

	~IntegratorInterface() = default;

	TensorTree<T> Derivative(double t, const TensorTree<T>& Psi) {
		::Derivative(dPsi_, hRep_, t, Psi, H_, tree_, phase_);
		return dPsi_;
	}
private:
	const SOP<T>& H_;
	const Tree& tree_;
	T phase_;
	tHamiltonianRepresentation<T> hRep_;
	TensorTree<T> dPsi_;

};

#endif //INTEGRATORINTERFACE_H
