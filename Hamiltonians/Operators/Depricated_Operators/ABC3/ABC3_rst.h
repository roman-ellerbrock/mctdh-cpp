#pragma once
#include "FortranSOP.h"

class ABC3_rst :
	public FortranSOP
{
public:
	ABC3_rst(const mctdhBasis &basis);
	ABC3_rst(const mctdhBasis &basis, Vectord m);
	~ABC3_rst();

	void callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag) override;

//	Vectord IntToCart(const Vectord& q) override;

	void InitOperator() override ;
private:
	Vectord mass;
};

