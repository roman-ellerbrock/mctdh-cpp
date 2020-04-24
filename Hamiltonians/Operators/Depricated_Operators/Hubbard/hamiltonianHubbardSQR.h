#pragma once
#include "Hamiltonian.h"
#include "mctdhMatrices.h"
#include "oSQRSPO.h"
#include "Expectationvalue.h"

class hamiltonianHubbardSQR : public Hamiltonian
{
 public:
	hamiltonianHubbardSQR(double J = 1., double U = 0., double freq_ = 0.);
	~hamiltonianHubbardSQR();

	void Initialize(const mctdhBasis & basis);

	void additionalDatOut(const mctdhWavefunction& Psi,
	                      const mctdhMatrices& matrices,
	                      const mctdhBasis& basis,
	                      ostream& os);

	void Correlation(const mctdhWavefunction& Psi,
	                 const mctdhMatrices& matrices,
	                 const mctdhBasis& basis,
	                 ostream& os);
	
	void Occupation(const mctdhWavefunction& Psi,
	                const mctdhMatrices& matrices,
	                const mctdhBasis& basis,
	                ostream& os);
protected:
	double U;
	double J;
	double freq;
};


void HubbardDatOut(const mctdhWavefunction& Psi,
                   const mctdhMatrices& matrices,
                   const mctdhBasis& basis,
                   ostream& os);

void HubbardCorrelation(const mctdhWavefunction& Psi,
                        const mctdhMatrices& matrices,
                        const mctdhBasis& basis,
                        ostream& os);
	
void HubbardOccupation(const mctdhWavefunction& Psi,
                       const mctdhMatrices& matrices,
                       const mctdhBasis& basis,
                       ostream& os);
