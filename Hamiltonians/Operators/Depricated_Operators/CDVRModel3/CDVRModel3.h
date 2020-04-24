#pragma once
#include "Hamiltonian.h"

class CDVRModel3 :
	public Hamiltonian
{
public:
	CDVRModel3();
	~CDVRModel3();

	Vectord IntToCart(const Vectord& q) { return q; }

	void SpecialInitialize(const mctdhBasis & basis);

private:
	bool simple;
};

