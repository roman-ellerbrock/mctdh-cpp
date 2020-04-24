#include "CDVRModel3.h"

CDVRModel3::CDVRModel3()
	:simple(false)
{
}

CDVRModel3::~CDVRModel3()
{
}

void CDVRModel3::SpecialInitialize(const mctdhBasis& basis)
{
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x2 = &PrimitiveBasis::ApplyX2;
	constexpr double cm = 219474.6313705;
	constexpr double lambda = 2000. / cm;
	constexpr double omega = 4000. / cm;
	constexpr double c = 0.5*omega*omega;
	size_t f = 4;

	cout << "=== Initializing CDVRModel Hamiltonian ===" << endl;

	for (size_t i = 0; i < f; i++)
	{
		MultiParticleOperator M;
		M.push_back(kin, i);
		hamiltonian.push_back(M);
		coeff.push_back(1.0);
	}

	{
		MultiParticleOperator M;
		M.push_back(x2, f - 1);
		hamiltonian.push_back(M);
		coeff.push_back(c);
	}

	if (simple) {
		{
			MultiParticleOperator M;
			M.push_back(x, f - 2);
			M.push_back(x, f - 1);
			hamiltonian.push_back(M);
			coeff.push_back(lambda*lambda);
		}
		{
			MultiParticleOperator M;
			M.push_back(x, f - 1);
			M.push_back(x, 0);
			hamiltonian.push_back(M);
			coeff.push_back(lambda*lambda);
		}
		cout << "Initialized operator with coupling parts." << endl;
	}else{
		cout << "Initialized operator with NO coupling parts." << endl;
	}

	cout << "=== End initializing Operator ===" << endl;

}




