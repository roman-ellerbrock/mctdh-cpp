//
// Created by Roman Ellerbrock on 4/3/20.
//

#ifndef EXPLICITEDGEWAVEFUNCTION_H
#define EXPLICITEDGEWAVEFUNCTION_H
#include "Core/Wavefunction.h"
#include "TreeClasses/MatrixTreeFunctions.h"

class ExplicitEdgeWavefunction : public pair<Wavefunction, MatrixTreecd> {
public:
	ExplicitEdgeWavefunction(const Wavefunction& Psi, const Tree& tree, bool orthogonal);
	/**
	 * \brief All-normalized wavefunction representation with A\tilde's and B_inv's.
	 */

	const TensorTreecd& nodes()const { return first; }
	TensorTreecd& nodes() { return first; }

	const MatrixTreecd& edges()const { return second; }
	MatrixTreecd& edges() { return second; }

	TensorTreecd TopDownNormalized(const Tree& tree) const;

	TensorTreecd BottomUpNormalized(const Tree& tree) const;

};


bool IsWorking_bottomup(const ExplicitEdgeWavefunction& Psi, const Tree& tree, double eps);
bool IsWorking_topdown(const ExplicitEdgeWavefunction& Psi, const Tree& tree, double eps);
bool IsWorking(const ExplicitEdgeWavefunction& Psi, const Tree& tree, double eps = 1e-7);


#endif //EXPLICITEDGEWAVEFUNCTION_H
