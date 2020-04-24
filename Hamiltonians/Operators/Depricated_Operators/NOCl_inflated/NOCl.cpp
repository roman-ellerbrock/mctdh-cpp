#include "NOCl.h"


NOCl::NOCl()
{
}

NOCl::~NOCl()
{
}

void NOCl::Initialize(const mctdhBasis& basis)
{
	/*
	InitV(basis, "../data/NOCl/NOClcoefficients.txt");
	assert(hamiltonian.size() == 21);
	assert(coeff.size() == 21);

	InitKin();
	assert(hamiltonian.size() == 25);
	assert(coeff.size() == 25);
	*/

	InitKin();
	assert(hamiltonian.size() == 4);
	assert(coeff.size() == 4);
	/*
	*/
}


// Parts for the moments of inertia
Tensorcd recrsq2(const PrimitiveBasis& grid, const Tensorcd& Acoeffs)
{

	const Vectord& x = grid.GetX();
	Tensorcd xAcoeffs(Acoeffs);
	for (int n = 0; n < Acoeffs.Dim().getntensor(); n++)
	{
		for (int i = 0; i < Acoeffs.Dim().getdimpart(); i++)
		{
			xAcoeffs(i, n) /= x(i)*x(i);
		}
	}
	return xAcoeffs;
}

void NOCl::InitKin()
{
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;

	SPO recrsqa = recrsq2;

	// kinetic energy
	// r
	{
		MultiParticleOperator M;
		M.push_back(kin, 0);

		coeff.push_back(1.0);
		hamiltonian.push_back(M);
	}
	// R
	{
		MultiParticleOperator M;
		M.push_back(kin, 1);

		coeff.push_back(1.0);
		hamiltonian.push_back(M);
	}

	// Theta
	{
		MultiParticleOperator M;
		M.push_back(recrsqa, 0);
		M.push_back(kin, 2);

		coeff.push_back(1.0);
		hamiltonian.push_back(M);
	}
	{
		MultiParticleOperator M;
		M.push_back(recrsqa, 1);
		M.push_back(kin, 2);

		coeff.push_back(1.0);
		hamiltonian.push_back(M);
	}
}

void NOCl::InitV(const mctdhBasis& basis, string filename)
{
	ifstream is(filename);
	// Theta
	const PhysicalCoordinate& phys = basis.Phys(2);
	const PrimitiveBasis& grid = phys.PrimitiveGrid();
	// R
	const PhysicalCoordinate& physDiss = basis.Phys(1);
	const PrimitiveBasis& gridDiss = physDiss.PrimitiveGrid();


	SPO ApplyNOpotentiala = ApplyNOpotential;

	for (int i = 0; i <= 3; i++)
	{
		for (int j = 0; j <= 4; j++)
		{
			MultiParticleOperator M;

			// q_vib - r
			if (i > 0)
			{
				NOCloperatorVib* qvib_ptr = new NOCloperatorVib;
				qvib_ptr->Initialize(i);
				shared_ptr<BottomLayerSPO> qvib(qvib_ptr);
				M.push_back(qvib, 0);
			}

			// q_diss - R
			NOCloperatorDissoziation* qd_ptr = new NOCloperatorDissoziation;
			qd_ptr->Initialize(j, gridDiss);
			shared_ptr<BottomLayerSPO> qd(qd_ptr);
			M.push_back(qd, 1);

			// q_theta - theta
			NOCloperatorTheta* qtheta_ptr = new NOCloperatorTheta;
			qtheta_ptr->Initialize(is, 7, grid);
			shared_ptr<BottomLayerSPO> qtheta(qtheta_ptr);
			M.push_back(qtheta, 2);

			// save operator
			hamiltonian.push_back(M);
			coeff.push_back(1.);
		}
	}

	// NO Potential
	{
		MultiParticleOperator M;
		M.push_back(ApplyNOpotentiala, 0);
		hamiltonian.push_back(M);
		coeff.push_back(1.);
	}

}
