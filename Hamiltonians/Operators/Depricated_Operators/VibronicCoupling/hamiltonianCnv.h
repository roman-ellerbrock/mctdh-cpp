#pragma once
#include "Hamiltonian.h"
#include "spoSB.h"
#include "mctdhDatout.h"

class hamiltonC5v : public Hamiltonian
{
public:
	hamiltonC5v(vector<double> qpcoeff,
	            vector<double> qmcoeff,
	            vector<double> vcoeff);
	~hamiltonC5v();

	void Initialize(const mctdhBasis & basis); //abitrary E x e jt
	void InitializeRelativistic(const mctdhBasis& basis); //only c3v

protected:
	void addcouplingterm(const mctdhBasis& basis, size_t order,
	                     double coeff, int sign);
	void addoffsetterm(const mctdhBasis& basis, size_t halforder,
	                   double coeff);
	double over(size_t n, size_t k);
	
	vector<double> qpcoeff;
	vector<double> qmcoeff;
	vector<double> vcoeff;
};


void CnvDatOut(const mctdhWavefunction& Psi,
               const mctdhMatrices& matrices,
               const mctdhBasis& basis,
               ostream& os);

void CnvAngularMomentum(const mctdhWavefunction& Psi,
                        const mctdhMatrices& matrices,
                        const mctdhBasis& basis,
                        ostream& os);

void CnvSpinState(const mctdhWavefunction& Psi,
                  const mctdhMatrices& matrices,
                  const mctdhBasis& basis,
                  ostream& os);
