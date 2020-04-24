#pragma once
#include "Hamiltonian.h"

class CH3_meanfield :
	public SOP {
public:
	CH3_meanfield() = default;
    CH3_meanfield(const mctdhBasis& basis) {
        Initialize(basis);
    }
	~CH3_meanfield() = default;

	void Initialize(const mctdhBasis & basis);

	void InitCH3meanfield(const mctdhBasis & basis);

protected:
	double functionGmeanfield(double x, int mode, int part);
};
