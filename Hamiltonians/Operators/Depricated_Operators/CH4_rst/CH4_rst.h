#pragma once
#include "FortranSOP.h"
class CH4_rst :
	public FortranSOP
{
public:
	CH4_rst();
	~CH4_rst();

	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	Vectord IntToCart(const Vectord& q);

	void InitOperator();
};
