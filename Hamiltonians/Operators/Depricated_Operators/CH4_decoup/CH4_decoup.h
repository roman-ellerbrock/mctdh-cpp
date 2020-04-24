#pragma once
#include "FortranSOP.h"
class CH4_decoup :
	public FortranSOP
{
public:
	CH4_decoup();
	~CH4_decoup();

	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	void InitOperator();
	
	Vectord IntToCart(const Vectord& q);
};
