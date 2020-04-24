#pragma once
#include "Hamiltonian.h"

class hamiltonKvN :
	public Hamiltonian
{
public:
	hamiltonKvN();
	~hamiltonKvN();

	void Initialize(const mctdhBasis & basis);
};

