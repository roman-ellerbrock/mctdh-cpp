//
// Created by Roman Ellerbrock on 3/1/20.
//

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "TreeClasses/TensorTree.h"

typedef TensorTreecd Wavefunction;

void occupyCIS(TensorTreecd& Psi, const Tree& tree);

#endif //WAVEFUNCTION_H
