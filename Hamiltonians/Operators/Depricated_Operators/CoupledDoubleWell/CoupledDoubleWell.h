//
// Created by Roman on 3/6/2019.
//

#ifndef MCTDH_COUPLEDDOUBLEWELL_H
#define MCTDH_COUPLEDDOUBLEWELL_H
#pragma once

#include "Hamiltonian.h"
#include "FermionNumberBasis.h"
#include "Pauli.h"

class CoupledDoubleWell :
	public SOP
{
public:
	CoupledDoubleWell(const mctdhBasis& basis, size_t f, int dim, double h, double J);
	~CoupledDoubleWell() = default;

	void SpecialInitialize(const mctdhBasis & basis);

	void InitializeIsing(const mctdhBasis& basis, size_t f, int dim, double h = 0., double J = 1.);

private:
	bool simple;
};



#endif //MCTDH_COUPLEDDOUBLEWELL_H
