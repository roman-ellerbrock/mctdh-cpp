#pragma once
#include "FortranHamiltonian.h"
class OHCH4 :
	public FortranHamiltonian
{
public:
	OHCH4();
	~OHCH4();

	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	Vectord IntToCart(const Vectord& q);

	void InitOperator();
};
