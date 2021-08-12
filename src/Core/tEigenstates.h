//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef TEIGENSTATES_H
#define TEIGENSTATES_H
#include "Core/tHamiltonianRepresentation.h"
#include "Core/tIntegratorInterface.h"
#include "TreeClasses/TreeIO.h"
#include "Core/tCMFIntegrator.h"

template<typename T>
Vectord tpropagatorEnergies(const TensorTree<T>& Psi, const Tree& tree, double out);

template<typename T>
Vectord tEigenstate(TensorTree<T>& Psi, const SOP<T>& H, const Tree& tree);

template<typename T>
void tStatus(const Vectord& eigenvalues, const Vectord& propergatorev, const Matrix<T>& S, ostream& os);

template<typename T>
void tEigenstates(tIntegratorVariables<T>& ivar);

#endif //TEIGENSTATES_H
