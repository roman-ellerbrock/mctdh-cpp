#pragma once
#include "TreeOperators/SumOfProductsOperator.h"

class Hamiltonian : public SOPcd
{
public:
	Hamiltonian() = default;
	~Hamiltonian() = default;

	using SOPcd::operator=;

};

