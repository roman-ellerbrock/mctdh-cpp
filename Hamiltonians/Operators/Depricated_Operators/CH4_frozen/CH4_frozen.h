#pragma once
#include "FortranSOP.h"
class CH4_frozen :
	public FortranSOP
{
public:
	CH4_frozen();
	~CH4_frozen();

	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	void InitOperator();
	
	Vectord IntToCart(const Vectord& q);
};
