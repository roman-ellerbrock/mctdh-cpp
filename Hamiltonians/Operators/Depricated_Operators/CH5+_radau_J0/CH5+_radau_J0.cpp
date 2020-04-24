#include "CH5+_radau_J0.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	void __stdcall HINITCH5P(double* coeff, int* diag, int* nmodes);
	// Hamiltonian in Fortran (SystemH)
	void __stdcall HCH5P(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
#ifdef __linux__
	void hinitch5p_(double* coeff, int* diag, int* nmodes);
	void hch5p_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
}


KEO_CH5Plus_J0::KEO_CH5Plus_J0()
{
}


KEO_CH5Plus_J0::~KEO_CH5Plus_J0()
{
}

void KEO_CH5Plus_J0::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	// Call Hinit here
#ifdef _MSC_VER
	HINITCH5P((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
#endif
#ifdef __linux__
	hinitch5p_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
#endif
}

void KEO_CH5Plus_J0::InitOperator()
{
	// set number of parts
	// k zahl=4+34*(N-3)+38*(N-3)*(N-4)+15*(N-3)+1+4
	// for N=6: 384
	// for J=0: k zahl=k zahl-15*(N-3)-4-1
	// which is 50 (for N=6)
	// -> (J=0, N=6) nparts=334, but need 384 for hinit
	nparts = 384;
	// number of physical modes
	nmodes = 12;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HCH5P;
#endif
#ifdef __linux__
	SystemH = hch5p_;
#endif
}


Vectord KEO_CH5Plus_J0::IntToCart(const Vectord& q)
{
	return q;
}
