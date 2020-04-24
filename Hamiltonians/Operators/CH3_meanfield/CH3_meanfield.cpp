#include "CH3_meanfield.h"

namespace ch3operators
{

//	constexpr double rhoeq = sqrt(3.*1836.109)*2.04035;
	constexpr double rhoeq = 74.218104260348768355888738518379*2.04035;
	constexpr double dphr2 = 0.0030;
	constexpr double dthr2 = 0.0045;
	constexpr double dph2 = 0.010;
	constexpr double dch2 = 0.030;  

	constexpr double pi = 3.1415926535897931;

	void ApplyV(Tensorcd& HPsi, const Tensorcd& Psi, const Vectord& x, function<double(double)> g, function<double(double)> f)
	{
		assert(HPsi.Dim() == Psi.Dim());
		assert(HPsi.Dim().getdimpart() == x.Dim());
		for (int n=0; n<Psi.Dim().getntensor(); n++)
			for (int i = 0; i < Psi.Dim().getdimpart(); i++)
				HPsi(i, n) = Psi(i, n)*g(f(x(i)));
	}

	double v_sqrt(double x)
	{
		assert(x > 0);
		return sqrt(x);
	}

	double v_divsqrt(double x)
	{
		assert(x > 0);
		return 1./sqrt(x);
	}

	void ApplyG(Tensorcd& HPsi, const Tensorcd& Psi, const Vectord& x, function<double(double)> f)
	{
		assert(HPsi.Dim() == Psi.Dim());
		assert(HPsi.Dim().getdimpart() == x.Dim());
		for (int n = 0; n < Psi.Dim().getntensor(); n++)
			for (int i = 0; i < Psi.Dim().getdimpart(); i++)
				HPsi(i, n) = Psi(i, n)*f(x(i));
	}

	void divx2(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) { return 1./(x*x); };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void a_sin(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return sin(x); };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void a_sin3(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return pow(sin(x), 3); };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void a_divsin2(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return 1./(sin(x)*sin(x)); };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	// V - rho
	void vol_rho(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return x*x; };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void divvol_rho(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) { return x*x; };
		ApplyV(HPsi, Psi, grid.GetX(), v_divsqrt, f);
	}

	// V - theta_rho
	void vol_thetarho(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		a_sin(grid, HPsi, Psi);
	}

	void divvol_thetarho(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return sin(x); };
		ApplyV(HPsi, Psi, grid.GetX(), v_divsqrt, f);
	}

	// V - theta
	void vol_theta(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		a_sin3(grid, HPsi, Psi);
	}

	void sqrtvol_theta(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return pow(sin(x), 3); };
		ApplyV(HPsi, Psi, grid.GetX(), v_sqrt, f);
	}

	void divvol_theta(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [](double x) {return pow(sin(x), 3); };
		ApplyV(HPsi, Psi, grid.GetX(), v_divsqrt, f);
	}

	// V - phi
	void vol_phi(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		constexpr double pi3 = pi / 3.;
		auto f = [=](double x)
		{
			return 1. / (1. + pow(x - pi3, 2) - pow(x - pi3, 3) / (3.*sqrt(3.)) + pow(x - pi3, 4)*3. / 4.);
		};
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void sqrtvol_phi(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		constexpr double pi3 = pi / 3.;
		auto f = [=](double x)
		{
			return 1. / (1. + pow(x - pi3, 2) - pow(x - pi3, 3) / (3.*sqrt(3.)) + pow(x - pi3, 4)*3. / 4.);
		};
		ApplyV(HPsi, Psi, grid.GetX(), v_sqrt, f);
	}

	void divvol_phi(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		constexpr double pi3 = pi / 3.;
		auto f = [=](double x)
		{
			return 1. / (1. + pow(x - pi3, 2) - pow(x - pi3, 3) / (3.*sqrt(3.)) + pow(x - pi3, 4)*3. / 4.);
		};
		ApplyV(HPsi, Psi, grid.GetX(), v_divsqrt, f);
	}

	// V - chi
	void vol_chi(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [=](double x)
		{
			double fun = 1. / (1. + pow(x - pi, 2) / 3. + pow(x - pi, 4) / 12.);
			return fun;
		};
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void divvol_chi(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [=](double x)
		{
			double fun = 1. / (1. + pow(x - pi, 2) / 3. + pow(x - pi, 4) / 12.);
			return fun;
		};
		ApplyV(HPsi, Psi, grid.GetX(), v_divsqrt, f);
	}

	// G
	void G91(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [=](double x) {return (1. + 2.*dphr2 + 4. / 3.*dthr2 + dph2 / 3. + dch2 / 9.); };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void G92(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [=](double x) {return 1. / (pow(tan(x), 2) * 4. * sqrt(3.)) *
			(-9.*dphr2 + 8.*dthr2 + 7.*dph2 - 7. / 3.*dch2); };
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void G93(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [=](double x)
		{
			return (3. / 4.*(3. + cos(2.*x)) / pow(sin(x), 2)
			+ 3. / 8.*(-11. + 16. / pow(sin(x), 2))*dphr2
			+ (-1. + 4. / pow(sin(x), 2))*dthr2
			+ 3. / 2. / pow(tan(x), 2) * dph2
			+ 5. / 6. / pow(tan(x), 2) * dch2);
		};
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}

	void G94(const PrimitiveBasis& grid, Tensorcd& HPsi, const Tensorcd& Psi)
	{
		auto f = [=](double x)
		{
			double g =
				( 9. / 4.*(3. + cos(2.*x)) / (sin(x)*sin(x))
				+ 9. / 8.*(-5. + 16. / pow(sin(x), 2))*dphr2
				+ 3. * (1. + 2. / pow(tan(x), 2))*dthr2
				+ 15. / 2. / pow(tan(x), 2) * dph2
				+ 3. / 2. / pow(tan(x), 2) * dch2);
			return g;
		};
		ApplyG(HPsi, Psi, grid.GetX(), f);
	}
}

void CH3_meanfield::Initialize(const mctdhBasis& basis) {
	InitCH3meanfield(basis);
}

void CH3_meanfield::InitCH3meanfield(const mctdhBasis& basis)
{
	using namespace ch3operators;

	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> kin = &PrimitiveBasis::ApplyKin;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> p = &PrimitiveBasis::ApplyP;

	// 0.5/sqrt(v)*p*v*p/sqrt(v)
	{
		MultiParticleOperator M;
		M.push_back(divvol_rho, 0);
		M.push_back(p, 0);
		M.push_back(vol_rho, 0);
		M.push_back(p, 0);
		M.push_back(divvol_rho, 0);
		
		push_back(M, 0.5);
	}

	// 1/rho�/sqrt(v)*p*v*p/sqrt(v)
	{
		MultiParticleOperator M;
		M.push_back(divx2, 0);

		M.push_back(divvol_thetarho, 1);
		M.push_back(p, 1);
		M.push_back(vol_thetarho, 1);
		M.push_back(p, 1);
		M.push_back(divvol_thetarho, 1);
		
		push_back(M, 0.5);
	}

	// part 3: 1/(rho�*sin�(phi_rho))*kin_theta_rho
	{
		MultiParticleOperator M;
		M.push_back(divx2, 0);
		M.push_back(a_divsin2, 1);

		// p�/2
		M.push_back(kin, 2);

		push_back(M, 1.);
	}

	// part 4
	// 1/rho�/sqrt(v)*p*G91*v_theta*p/sqrt(v)
	{
		MultiParticleOperator M;
		M.push_back(divx2, 0);

		M.push_back(divvol_theta, 3);
		M.push_back(p, 3);
		M.push_back(G91, 3);
		M.push_back(vol_theta, 3);
		M.push_back(p, 3);
		M.push_back(divvol_theta, 3);
		
		push_back(M, 0.5);
	}

	// part 5
	{
		MultiParticleOperator M;
		M.push_back(divx2, 0);

		M.push_back(G92, 3);
		M.push_back(sqrtvol_theta, 3);
		M.push_back(p, 3);
		M.push_back(divvol_theta, 3);

		M.push_back(divvol_phi, 4);
		M.push_back(p, 4);
		M.push_back(sqrtvol_phi, 4);

		push_back(M, 0.5);

	}

	// Part 6
	{
		// 1/rho�
		MultiParticleOperator M;
		M.push_back(divx2, 0);

		// Theta
		M.push_back(divvol_theta, 3);
		M.push_back(p, 3);
		M.push_back(sqrtvol_theta, 3);
		M.push_back(G92, 3);

		// Phi
		M.push_back(sqrtvol_phi, 4);
		M.push_back(p, 4);
		M.push_back(divvol_phi, 4);

		push_back(M, 0.5);
	}

	// Part 7
	{
		MultiParticleOperator M;
		// 1/rho�
		M.push_back(divx2, 0);

		// G in theta
		M.push_back(G93, 3);

		// Phi
		M.push_back(divvol_phi, 4);
		M.push_back(p, 4);
		M.push_back(vol_phi, 4);
		M.push_back(p, 4);
		M.push_back(divvol_phi, 4);

		push_back(M, 0.5);
	}

	// Part 8
	{
		MultiParticleOperator M;
		// 1/rho�
		M.push_back(divx2, 0);

		// G in theta
		M.push_back(G94, 3);

		// Chi
		M.push_back(divvol_chi, 5);
		M.push_back(p, 5);
		M.push_back(vol_chi, 5);
		M.push_back(p, 5);
		M.push_back(divvol_chi, 5);

		push_back(M, 0.5);
	}

}
