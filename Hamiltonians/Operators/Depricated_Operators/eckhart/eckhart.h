#pragma once
#include "FortranHamiltonian.h"
class eckhart :
	public FortranHamiltonian
{
public:
	eckhart();
	~eckhart();

	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	void InitOperator();
};
