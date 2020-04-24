#include "eckhart.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	// Hamiltonian in Fortran (SystemH)
	void __stdcall HINITECKHART(double* coeff, int* diag, int* nmodes);
	void __stdcall HECKHART(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
#ifdef __linux__
	void hiniteckhart_(double* coeff, int* diag, int* nmodes);
	void heckhart_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
}

eckhart::eckhart()
{
}

eckhart::~eckhart()
{
}

void eckhart::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	// Call Hinit here
#ifdef _MSC_VER
	HINITeckhart((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
#endif
#ifdef __linux__
	hiniteckhart_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
#endif
}

void eckhart::InitOperator()
{
	// set number of parts
	nparts = 4;
	// number of physical modes
	nmodes = 2;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HECKHART;
#endif
#ifdef __linux__
	SystemH = heckhart_;
#endif
}

