#pragma once
#include "FortranSOP.h"
class HCH4_rst :
	public FortranSOP
{
public:
	HCH4_rst();
	~HCH4_rst();

	void callHinit(Vectorcd& coeffs, Matrix<int>& diag);
	void InitOperator();
	Vectord IntToCart(const Vectord& q);
};

