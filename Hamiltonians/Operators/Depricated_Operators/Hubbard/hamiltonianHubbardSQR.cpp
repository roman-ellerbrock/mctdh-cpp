#include "hamiltonianHubbardSQR.h"

hamiltonianHubbardSQR::hamiltonianHubbardSQR(double J_, double U_, 
                                             double freq_)
{
	J = J_;
	U = U_;
	freq = freq_;
}

hamiltonianHubbardSQR::~hamiltonianHubbardSQR(){}

void hamiltonianHubbardSQR::Initialize(const mctdhBasis& basis)
{	
	hamiltonian.clear();
	coeff.clear();

	//Parameters for particle number conservation
	double strength = 10000;
	double decay = 0.1;

	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		x = &PrimitiveBasis::applyX;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		p = &PrimitiveBasis::ApplyP;

	//get number of modes
	size_t prim = basis.nPhysNodes();
	
	//Kinetic energy
	for(size_t i = 0; i < prim; i++)
	{
		{
			MultiParticleOperator M;
		
			M.push_back(x,i);
			if(i+1 < prim)
			{
				M.push_back(p,i+1);
			}
			else
			{
				M.push_back(p,i+1-prim);
			}
			
			coeff.push_back(-J);
			hamiltonian.push_back(M);
		}
		{
			MultiParticleOperator M;
		
			if(i+1 < prim)
			{
				M.push_back(x,i+1);
			}
			else
			{
				M.push_back(x,i+1-prim);
			}
			M.push_back(p,i);
		
			coeff.push_back(-J);
			hamiltonian.push_back(M);
		}
	}

	//On side repulsion
	for(size_t i = 0; i < prim; i++)
	{
		MultiParticleOperator M;
		
		M.push_back(x,i);
		M.push_back(x,i);
		M.push_back(p,i);
		M.push_back(p,i);
		
		coeff.push_back(U);
		hamiltonian.push_back(M);
	}

	//harmonic trap
	for(size_t i = 0; i < prim; i++)
	{
		MultiParticleOperator M;
		
		M.push_back(x,i);
		M.push_back(p,i);
		
		coeff.push_back(0.5*pow(freq*(1.*i-0.5*(prim-1)), 2));
		hamiltonian.push_back(M);
	}

	//Particle Conservation
	{
		MultiParticleOperator M1;
		MultiParticleOperator M2;
		MultiParticleOperator M3;
		
		//get particle number
		size_t N = 0;
		for(int i = 0; i < prim; i++)
		{
			const PhysicalCoordinate& phys = basis.Phys(i);
			
			PhysPar par = phys.Par();
			int occupation = par.Omega();
		  
			N += occupation;
		}
		
		for(int i = 0; i < prim; i++)
		{
			const PhysicalCoordinate& phys = basis.Phys(i);
		
			PhysPar par = phys.Par();
			int occupation = par.Omega();
			int minOcc = par.WFR0();
	  
			{
				shared_ptr<BottomLayerSPO> conserver
					= shared_ptr<BottomLayerSPO>(new oSQRSPO(decay/sqrt(N),
					                                         minOcc));
	
				M1.push_back(move(conserver), i);
			}
				
			{
				shared_ptr<BottomLayerSPO> conserver
					= shared_ptr<BottomLayerSPO>(new oSQRSPO(-decay/sqrt(N),
					                                         minOcc));
				
				M2.push_back(move(conserver), i);
			}
		}
		
		coeff.push_back(strength*N*exp(-decay*sqrt(N)));
		coeff.push_back(strength*N*exp(decay*sqrt(N)));
		coeff.push_back(-2.*strength*N);
	
		hamiltonian.push_back(move(M1));
		hamiltonian.push_back(move(M2));
		hamiltonian.push_back(move(M3));
	}

	//Zeropoint Energy
	{
		//get particle number
		size_t N = 0;
		for(int i = 0; i < prim; i++)
		{
			const PhysicalCoordinate& phys = basis.Phys(i);
			
			PhysPar par = phys.Par();
			int occupation = par.Omega();
		  
			N += occupation;
		}

		{
			MultiParticleOperator M;
		
			coeff.push_back(2.*J*N);
			hamiltonian.push_back(M);
		}
	}
}

void HubbardDatOut(const mctdhWavefunction& Psi,
                   const mctdhMatrices& matrices,
                   const mctdhBasis& basis,
                   ostream& os)
{
	HubbardCorrelation(Psi, matrices, basis, os);
	HubbardOccupation(Psi, matrices, basis, os);
}

void HubbardCorrelation(const mctdhWavefunction& Psi,
                        const mctdhMatrices& matrices,
                        const mctdhBasis& basis,
                        ostream& os)
{
	using namespace expectationvalue;

	//get ladder operators (fermionic or bosonic)
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		kreator = &PrimitiveBasis::ApplyP;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		annihilator = &PrimitiveBasis::applyX;
		
	//for results
	vector< complex<double> > corrMat;

	//Number of coordinates
	double act = basis.nPhysNodes();

	//calculate all two point functions
	for(size_t i = 0; i < act; i++)
	{
		for(size_t j = 0; j < act; j++)
		{
			//make a^+_j a_i
			MultiParticleOperator M;
			M.push_back(annihilator,i);
			M.push_back(kreator,j);

			corrMat.push_back(StateExpectationvalue(M,Psi,basis,0));
		}
	}

	//print the values of the two point functions in gnuplot input format
	os << endl;
	os << "two point functions" << endl;
	for(size_t i = 0; i < act; i++)
	{
		for(size_t j = 0; j < act; j++)
		{
			os << i << " " << j << " " << real(corrMat[i*act+j])
			   << " " << imag(corrMat[i*act+j]) << endl;
		}
		os << endl;
	}	
}

void HubbardOccupation(const mctdhWavefunction& Psi,
                       const mctdhMatrices& matrices,
                       const mctdhBasis& basis,
                       ostream& os)
{
	using namespace expectationvalue;

	//get ladder operators (fermionic or bosonic)
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		kreator = &PrimitiveBasis::ApplyP;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)>
		annihilator = &PrimitiveBasis::applyX;
		
	//for results
	vector< complex<double> > occupation;

	//Number of coordinates
	double act = basis.nPhysNodes();
	
	//calculate all occupations
	for(size_t i = 0; i < act; i++)
	{
		//make a^+_j a_i
		MultiParticleOperator M;
		M.push_back(annihilator,i);
		M.push_back(kreator,i);

		occupation.push_back(StateExpectationvalue(M,Psi,basis,0));
	}

	//print the occupations in gnuplot input format
	os << endl;
	os << "Occupations" << endl;
	for(size_t i = 0; i < act; i++)
	{
		os << i << " " << abs(occupation[i]) << endl;
	}	
	os << endl;
}
