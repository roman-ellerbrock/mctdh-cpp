#pragma once
#include "Potential.h"


class PES_CH3Cl :
	public Potential
{
public:
	PES_CH3Cl() = default;
	~PES_CH3Cl() = default;

	vector<MultiParticleOperator> Build(const mctdhBasis& basis);

	double Evaluate(const Vectord& Xv, size_t part);

protected:
	

};

