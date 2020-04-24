#include "hamiltonKvN.h"

hamiltonKvN::hamiltonKvN(){}

hamiltonKvN::~hamiltonKvN(){}

void hamiltonKvN::Initialize(const mctdhBasis& basis)
{
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;

	{
		MultiParticleOperator M;
		
		M.push_back(x,1);
		M.push_back(p,0);

		coeff.push_back(-1.);
		hamiltonian.push_back(M);
	}

	{
		MultiParticleOperator M;
		
		M.push_back(x,0);
		M.push_back(p,1);

		coeff.push_back(-1.);
		hamiltonian.push_back(M);
	}

	{
		MultiParticleOperator M;
		
		M.push_back(x,0);
		M.push_back(x,0);
		M.push_back(x,0);
		M.push_back(p,1);

		coeff.push_back(0.4);
		hamiltonian.push_back(M);
	}
}
