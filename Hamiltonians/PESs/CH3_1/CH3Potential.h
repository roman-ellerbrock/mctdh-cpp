#pragma once
#include "TreeOperators/Potential.h"

class CH3Potential :
	public Potential
{
public:
	CH3Potential() = default;
	~CH3Potential() = default;

	double Evaluate(const Vectord& Xv, size_t part) const override;
};

