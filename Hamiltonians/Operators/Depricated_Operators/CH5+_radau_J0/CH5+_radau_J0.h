#pragma once
#include "FortranSOP.h"
class KEO_CH5Plus_J0 :
	public FortranSOP
{
public:
	KEO_CH5Plus_J0();
	~KEO_CH5Plus_J0();

	// Call Fortran H Init
	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	// Set Variables for interface
	void InitOperator();

	// Internal coordinates To Cartesians
	Vectord IntToCart(const Vectord& q);
};

