//
// Created by Roman Ellerbrock on 4/3/20.
//

#ifndef EXPLICITEDGEWAVEFUNCTION_H
#define EXPLICITEDGEWAVEFUNCTION_H
#include "Core/Wavefunction.h"
#include "TreeClasses/MatrixTreeFunctions.h"

class ExplicitEdgeWavefunction : public pair<Wavefunction, MatrixTreecd> {
public:
	ExplicitEdgeWavefunction(const Wavefunction& Psi, const Tree& tree);
};


#endif //EXPLICITEDGEWAVEFUNCTION_H
