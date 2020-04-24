#pragma once
#include "FortranSOP.h"
class KEO_CH4Rad :
	public FortranSOP
{
public:
	KEO_CH4Rad();
	~KEO_CH4Rad();

	// Call Fortran H Init
	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag);

	// Coordinate transformation
	Vectord IntToCart(const Vectord& q);

	// Set Variables for interface
	void InitOperator();
};

