#include "hamiltonianCnv.h"

Tensorcd zerofield(const PrimitiveBasis& grid,const Tensorcd& phi)
{
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer tensor
	assert(tdim.F() == 1);

	Tensorcd psi(phi.Dim());

	complex<double> im(0.,1.);
	for (int n = 0; n < tdim.getntensor(); n++)
	{
		psi(1, n) = im*phi(0, n);
		psi(0, n) = -im*phi(1, n);
		psi(3, n) = im*phi(2, n);
		psi(2, n) = -im*phi(3, n);
	}
	return psi;
}

Tensorcd sorecoupling(const PrimitiveBasis& grid,const Tensorcd& phi)
{
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer tensor
	assert(tdim.F() == 1);

	Tensorcd psi(phi.Dim());

	for (int n = 0; n < tdim.getntensor(); n++)
	{
		psi(0, n) = phi(2, n);
		psi(1, n) = -phi(3, n);
		psi(2, n) = phi(0, n);
		psi(3, n) = -phi(1, n);
	}
	return psi;
}

Tensorcd soimcoupling(const PrimitiveBasis& grid,const Tensorcd& phi)
{
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer tensor
	assert(tdim.F() == 1);

	Tensorcd psi(phi.Dim());

	complex<double> im(0.,1.);
	for (int n = 0; n < tdim.getntensor(); n++)
	{
		psi(0, n) = -im*phi(2, n);
		psi(1, n) = im*phi(3, n);
		psi(2, n) = im*phi(0, n);
		psi(3, n) = -im*phi(1, n);
	}
	return psi;
}

Tensorcd jtwcoupling(const PrimitiveBasis& grid,const Tensorcd& phi)
{
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer tensor
	assert(tdim.F() == 1);

	Tensorcd psi(phi.Dim());

	for (int n = 0; n < tdim.getntensor(); n++)
	{
		psi(0, n) = -phi(0, n);
		psi(1, n) = phi(1, n);
		psi(2, n) = phi(2, n);
		psi(3, n) = -phi(3, n);
	}
	return psi;
}

Tensorcd jtzcoupling(const PrimitiveBasis& grid,const Tensorcd& phi)
{
	TensorDim tdim = phi.Dim();
	// check that its really a bottom-layer tensor
	assert(tdim.F() == 1);

	Tensorcd psi(phi.Dim());

	for (int n = 0; n < tdim.getntensor(); n++)
	{
		psi(0, n) = phi(1, n);
		psi(1, n) = phi(0, n);
		psi(2, n) = phi(3, n);
		psi(3, n) = phi(2, n);
	}
	return psi;
}

hamiltonC5v::hamiltonC5v(vector<double> qpcoeff_,
                         vector<double> qmcoeff_,
                         vector<double> vcoeff_)
{
	qpcoeff = qpcoeff_;
	qmcoeff = qmcoeff_;
	vcoeff = vcoeff_;
}

hamiltonC5v::~hamiltonC5v(){}

using namespace sposb;
void hamiltonC5v::Initialize(const mctdhBasis& basis)
{
	//Units: h/2pi = 1, m = 1
	//Koordiantes: x, y, elektronic
	
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		p = &PrimitiveBasis::ApplyP;

	SPO sigmax = &sposb::pauliX;
	SPO sigmaz = &sposb::pauliZ;

	//Energy offset
	{
		{
			MultiParticleOperator M;
	  	hamiltonian.push_back(M);
			coeff.push_back(35.);
		}
	}

	//Kinetic energy
	{
		{
			MultiParticleOperator M;
		
			M.push_back(p,0);
			M.push_back(p,0);

			coeff.push_back(0.5);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
		
			M.push_back(p,1);
			M.push_back(p,1);

			coeff.push_back(0.5);
			hamiltonian.push_back(M);
		}
	}

	//potential
	for(size_t i = 0; i < vcoeff.size(); i++)
	{
		double co = vcoeff[i];
		if(co != 0.)
		{
			addoffsetterm(basis, i+1, co);
		}
	}

	//field splitting
	for(size_t i = 0; i < qpcoeff.size(); i++)
	{
		double co = qpcoeff[i];
		if(co != 0.)
		{
			addcouplingterm(basis, i+1, co,  1);
		}
	}
	for(size_t i = 0; i < qmcoeff.size(); i++)
	{
		double co = qmcoeff[i];
		if(co != 0.)
		{
			addcouplingterm(basis, i+1, co, -1);
		}
	}
}

void hamiltonC5v::addcouplingterm(const mctdhBasis& basis, size_t order,
                                  double co, int sign)
{
	assert(abs(sign) == 1);
	assert(order > 0);
	
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		p = &PrimitiveBasis::ApplyP;

	SPO sigmax = &sposb::pauliX;
	SPO sigmaz = &sposb::pauliZ;

	//Diagonal elements
	{
		double rel = -1.;
		for(int i = 0; i <= order; i += 2)
		{
			MultiParticleOperator M;

			for(int j = 0; j < order-i; j++)
				M.push_back(x,0);			
			for(int j = 0; j < i; j++)
				M.push_back(x,1);
			
			M.push_back(sigmaz,2);
			
			double binomial = over(order, i);
			rel *= -1.;
			
			coeff.push_back(co*binomial*rel);
			hamiltonian.push_back(M);
		}
	}

	//off diagonal elements
	{
		double rel = -1.*sign;
		for(int i = 1; i <= order; i += 2)
		{
			MultiParticleOperator M;
			
			for(int j = 0; j < order-i; j++)
				M.push_back(x,0);
			for(int j = 0; j < i; j++)
				M.push_back(x,1);
			
			M.push_back(sigmax,2);
			
			double binomial = over(order, i);
			rel *= -1.;
			
			coeff.push_back(co*binomial*rel);
			hamiltonian.push_back(M);
		}
	}
}

void hamiltonC5v::addoffsetterm(const mctdhBasis& basis, size_t order,
                                double co)
{
	assert(order > 0);
	
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		p = &PrimitiveBasis::ApplyP;

	for(int i = 0; i <= order; i++)
	{
		MultiParticleOperator M;

		for(int j = 0; j < 2*(order-i); j++)
			M.push_back(x,0);			
		for(int j = 0; j < 2*i; j++)
			M.push_back(x,1);
			
		double binomial = over(order, i);
			
		coeff.push_back(co*binomial);
		hamiltonian.push_back(M);
	}
}


double hamiltonC5v::over(size_t n, size_t k)
{
	if(k > n) return 0.;
	size_t l = n-k;
	if(l > k)
	{
		size_t i;
		i = k;
		k = l;
		l = i;
	}

	double result = 1.;
	for(size_t i = n; i > k; i--)
	{
		result *= 1.*i;
		while(l > 0 && result > l)
		{
			result /= 1.*l;
			l--;
		}
	cout << result << endl;
	}
	while(l > 0)
	{
		result /= 1.*l;
		l--;
	}
	
	return result;
}

void hamiltonC5v::InitializeRelativistic(const mctdhBasis& basis)
{
	//Units: h/2pi = 1, m = 1
	//Koordiantes: x, y, elektronic

	double a = 3.0;//0.4; //first order jt coupling
	double epsilon = 0.5;//1; //Zero field splitting
	double theta = 2.;//0.3; //first order so coupling
	
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		p = &PrimitiveBasis::ApplyP;
	
	//Kinetic energy
	{
		{
			MultiParticleOperator M;
		
			M.push_back(p,0);
			M.push_back(p,0);

			coeff.push_back(0.5);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
		
			M.push_back(p,1);
			M.push_back(p,1);

			coeff.push_back(0.5);
			hamiltonian.push_back(M);
		}
	}
	
	//potential
	{
		{
			MultiParticleOperator M;
		
			M.push_back(x,0);
			M.push_back(x,0);

			coeff.push_back(0.5);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
		
			M.push_back(x,1);
			M.push_back(x,1);

			coeff.push_back(0.5);
			hamiltonian.push_back(M);
		}
	}

	//diabatic model
	{
		{
			MultiParticleOperator M;
			SPO zf = zerofield;

			M.push_back(zf,2);

			coeff.push_back(epsilon);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
			SPO sore = sorecoupling;

			M.push_back(x,1);
			M.push_back(sore,2);

			coeff.push_back(theta);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
			SPO soim = soimcoupling;

			M.push_back(x,0);
			M.push_back(soim,2);

			coeff.push_back(theta);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
			SPO jtw = jtwcoupling;

			M.push_back(x,0);
			M.push_back(jtw,2);

			coeff.push_back(a);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
			SPO jtz = jtzcoupling;

			M.push_back(x,1);
			M.push_back(jtz,2);

			coeff.push_back(a);
			hamiltonian.push_back(M);
		}
	}
}

void CnvDatOut(const mctdhWavefunction& Psi,
               const mctdhMatrices& matrices,
               const mctdhBasis& basis,
               ostream& os)
{
	CnvAngularMomentum(Psi, matrices, basis, os);
	CnvSpinState(Psi, matrices, basis, os);

	mctdhDatout Datout;
	Datout.Dump2D(Psi, basis, 0, 1);
}

void CnvAngularMomentum(const mctdhWavefunction& Psi,
                        const mctdhMatrices& matrices,
                        const mctdhBasis& basis,
                        ostream& os)
{
	using namespace expectationvalue;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		p = &PrimitiveBasis::ApplyP;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		x = &PrimitiveBasis::applyX;

	vector< complex<double> > angularMomentum;
	
	MultiParticleOperator M;
	M.push_back(p,1);
	M.push_back(x,0);
	MultiParticleOperator N;
	M.push_back(p,0);
	M.push_back(x,1);

	const mctdhNode& node = basis.TopNode();
	const TensorDim& dim = node.TDim();
	size_t states = dim.getntensor();

	for(size_t i = 0; i < states; i++)
	{
		complex<double> expect = StateExpectationvalue(M, Psi, basis, i);
		expect += StateExpectationvalue(N, Psi, basis, i);
		angularMomentum.push_back(expect);
	}
	
	cout << endl;
	cout << "Angular Momentum of spacial coordinates";
	for(size_t i = 0; i < states; i++)
		cout << angularMomentum[i];
	cout << endl;
	cout << endl;
}

void CnvSpinState(const mctdhWavefunction& Psi,
                  const mctdhMatrices& matrices,
                  const mctdhBasis& basis,
                  ostream& os)
{
	using namespace expectationvalue;
	const mctdhNode& node = basis.TopNode();
	const TensorDim& dim = node.TDim();
	size_t states = dim.getntensor();
	
	vector< complex<double> > angularMomentum;

	using namespace sposb;
	SPO sigmaz = &sposb::pauliZ;

	MultiParticleOperator S;
	S.push_back(sigmaz,2);

	for(size_t i = 0; i < states; i++)
		angularMomentum.push_back(StateExpectationvalue(S, Psi, basis, i));

	cout << endl;
	cout << "Angular Momentum of electronic state (or spin)";
	for(size_t i = 0; i < states; i++)
		cout << angularMomentum[i];
	cout << endl;
	cout << endl;
}
