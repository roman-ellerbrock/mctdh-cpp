//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef EIGENSTATES_H
#define EIGENSTATES_H
#include "Core/HamiltonianRepresentation.h"
#include "Core/IntegratorInterface.h"
#include "TreeClasses/TreeIO.h"
#include "Core/CMFIntegrator.h"

Vectord propagatorEnergies(const Wavefunction& Psi, const Tree& tree, double out);

Vectord Eigenstate(Wavefunction& Psi, const Hamiltonian& H, const Tree& tree);

void Status(const Vectord& eigenvalues, const Vectord& propergatorev, const Matrixcd& S, ostream& os);

void Eigenstates(IntegratorVariables& ivar);

#endif //EIGENSTATES_H
