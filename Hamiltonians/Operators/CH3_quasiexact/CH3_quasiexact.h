#pragma once
#include "FortranSOP.h"

class CH3_quasiexact :
	public FortranSOP
{
public:
	CH3_quasiexact(const mctdhBasis& basis);
	~CH3_quasiexact() = default;

private:
	void callHinit(Vectorcd& coeffs, Matrix<int>& diag) override;

	void InitOperator() override;
};

