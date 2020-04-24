#pragma once
#include "nocloperators.h"
#include "Hamiltonian.h"
#include "NOCl.h"
#include "NOClPotential.h"

class NOCl :
	public Hamiltonian
{
public:
	NOCl();
	~NOCl();

	void Initialize(const mctdhBasis & basis);

protected:
	void InitKin();
	void InitV(const mctdhBasis & basis, string filename);

};

