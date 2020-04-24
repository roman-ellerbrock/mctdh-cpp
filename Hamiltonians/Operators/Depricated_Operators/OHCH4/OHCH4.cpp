#include "OHCH4.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	// Hamiltonian in Fortran (SystemH)
	void __stdcall HINITOHCH4(double* coeff, int* diag, int* nmodes);
	// Initialization of masses for the coordinate transformation
	void __stdcall KOORDREADOHCH4();
	// Coordinate Transformation
	void __stdcall INTTOCARTOHCH4(double* q, double* x);
	void __stdcall HOHCH4(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
#ifdef __linux__
	void hinitohch4_(double* coeff, int* diag, int* nmodes);
	void hohch4_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void inttocartohch4_(double* q, double* x);
	void koordreadohch4_();
#endif
}

OHCH4::OHCH4()
{
}

OHCH4::~OHCH4()
{
}

void OHCH4::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	// Call Hinit here
#ifdef _MSC_VER
	HINITCH4RST((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	KOORDREADCH4RST();
#endif
#ifdef __linux__
	hinitch4rst_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	koordreadch4rst_();
#endif
}

void OHCH4::InitOperator()
{
	// set number of parts
	nparts = 69;
	// number of physical modes
	nmodes = 9;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HOHCH4;
#endif
#ifdef __linux__
	SystemH = hohch4_;
#endif
}

Vectord OHCH4::IntToCart(const Vectord& q)
{
	Vectord qq(q);
	Vectord Xv(15);
#ifdef _MSC_VER
	INTTOCARTOHCH4((double*) (&qq(0)), (double*) &Xv(0));
#endif
#ifdef __linux__
	inttocartohch4_((double*) (&qq(0)), (double*) &Xv(0));
#endif

	return Xv;
}
