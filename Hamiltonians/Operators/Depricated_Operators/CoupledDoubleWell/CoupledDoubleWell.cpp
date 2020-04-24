#include "CoupledDoubleWell.h"

Tensorcd ApplyWell(const PrimitiveBasis& grid,
		const Tensorcd& Acoeffs, int exponent)
{
	const Vectord& x = grid.GetX();
	Tensorcd xAcoeffs(Acoeffs);

	// mu stuff
	const double mh = 1836.;
	const double mN = 14.0*mh;
	const double mCl = 35.403*mh;
	const double mO = 16.0*mh;
	const double muv = 1. / mN + 1. / mO;
	const double mud = 1. / (mN + mO) + 1. / mCl;

	// equilibrium geometry
	const double alpha = 1.5;
	const double rde = 4.315;
	for (int i = 0; i < Acoeffs.Dim().getdimpart(); i++)
	{
		double factor = 1 - exp(-alpha*(x(i)*sqrt(mud) - rde));
		factor = pow(factor, exponent);
		for (int n = 0; n < Acoeffs.Dim().getntensor(); n++)
		{
			xAcoeffs(i, n) *= factor;
		}
	}
	return xAcoeffs;
}

CoupledDoubleWell::CoupledDoubleWell(const mctdhBasis& basis, size_t f, int dim, double h, double J)
	:simple(true) {
//	SpecialInitialize(basis);
	InitializeIsing(basis, f, dim, h, J);
}

void CoupledDoubleWell::SpecialInitialize(const mctdhBasis& basis)
{
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x2 = &PrimitiveBasis::ApplyX2;
	constexpr double cm = 219474.6313705;
	constexpr double a = 1.0/(16*1.3544);
	constexpr double b = -0.5;
	constexpr double c = 0.5;
	size_t f = 5;
	bool star = true;
	bool ring = false;

	cout << "=== Initializing Double-Well Hamiltonian ===" << endl;
	assert(0);

	for (size_t i = 0; i < f; i++)
	{
		MultiParticleOperator M;
		M.push_back(kin, i);
		mpos.push_back(M);
		coeff.push_back(1.0);
	}

	for (size_t i = 0; i < f; i++)
	{
		MultiParticleOperator M;
		M.push_back(x2, i);
		M.push_back(x2, i);
		mpos.push_back(M);
		coeff.push_back(a);
	}

	for (size_t i = 0; i < f; i++)
	{
		MultiParticleOperator M;
		M.push_back(x2, i);
		mpos.push_back(M);
		coeff.push_back(b);
	}

	if (ring) {
		for (size_t i = 0; i < f; i++) {
			MultiParticleOperator M;
			M.push_back(x, i);
			M.push_back(x, (i + 1)% f);
			mpos.push_back(M);
			coeff.push_back(c);
		}
	}

	if (star) {
		for (size_t i = 1; i < f; i++) {
			MultiParticleOperator M;
			M.push_back(x, 0);
			M.push_back(x, i);
			mpos.push_back(M);
			coeff.push_back(c);
		}
	}
	cout << "=== End initializing Operator ===" << endl;

}

void CoupledDoubleWell::InitializeIsing(const mctdhBasis& basis, size_t f, int dim, double h, double J) {
	using namespace PauliMatrices;

	cout << "Initializing Ising-operator with J=" << J << " and magnetic field of h=" << h << "\n";

	if (dim == 1) {
		cout << "Using " << f << " 1/2 spins on a 1D lattice.\n";
		if (h != 0.0) {
			for (size_t i = 0; i < f; ++i) {
				MultiParticleOperator M;
				M.push_back(sigma_x, i);
				mpos.push_back(M);
				coeff.push_back(-h);
			}
		}

		// 1D grid (periodic)
		for (size_t i = 0; i < f; ++i) {
			MultiParticleOperator M;
			M.push_back(sigma_z, i);
			M.push_back(sigma_z, (i + 1) % f);
			mpos.push_back(M);
			coeff.push_back(-J);
		}
	}else if (dim == 2) {
		cout << "2D lattice of 1/2 spins with tight binding interaction. Total number of spins is : " << f*f << endl;
		if (h != 0.0) {
			for (size_t i = 0; i < f; ++i) {
				for (size_t j = 0; j < f; ++j) {
				MultiParticleOperator M;
				M.push_back(sigma_x, idx(i, j, f));
				mpos.push_back(M);
				coeff.push_back(-h);
				}
			}
		}
		for (size_t i = 0; i < f; ++i) {
			for (size_t j = 0; j < f; ++j) {
				cout << i << " "  << j << " "  << idx(i, j, f) << endl;
				// top
				  {
				MultiParticleOperator M;
				M.push_back(sigma_z, idx(i, j, f));
				M.push_back(sigma_z, idx(i, j - 1, f));
				mpos.push_back(M);
				coeff.push_back(-J);
				  }

				// bottom
				  {
				MultiParticleOperator M;
				M.push_back(sigma_z, idx(i, j, f));
				M.push_back(sigma_z, idx(i, j + 1, f));
				mpos.push_back(M);
				coeff.push_back(-J);
				  }

				// right
				  {
				MultiParticleOperator M;
				M.push_back(sigma_z, idx(i, j, f));
				M.push_back(sigma_z, idx(i + 1, j, f));
				mpos.push_back(M);
				coeff.push_back(-J);
				  }
				
				// left
				  {
				MultiParticleOperator M;
				M.push_back(sigma_z, idx(i, j, f));
				M.push_back(sigma_z, idx(i - 1, j, f));
				mpos.push_back(M);
				coeff.push_back(-J);
				  }
			}
		}
	}else if (dim == -1) {
		cout << f << " 1/2 spins in infinit number of degrees of freedom (all spins coupled)" << endl;
		if (h != 0.0) {
			for (size_t i = 0; i < f; ++i) {
				MultiParticleOperator M;
				M.push_back(sigma_x, i);
				mpos.push_back(M);
				coeff.push_back(-h);
			}
		}
		for (size_t i = 0; i < f; ++i) {
			for (size_t j = i + 1; j < f; ++j) {
				MultiParticleOperator M;
				M.push_back(sigma_z, i);
				M.push_back(sigma_z, j);
				mpos.push_back(M);
				coeff.push_back(-J);
			}
		}
	}
}


