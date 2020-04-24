#include "HCH4_rst.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	void __stdcall HINITHCH4RST(double* coeff, int* diag, int* nmodes);
#endif
#ifdef __linux__
	void hinithch4rst_(double* coeff, int* diag, int* nmodes);
	void hhch4rst_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void inttocarthch4rst_(double* q, double* x);
	void koordreadhch4rst_();
#endif
	// Hamiltonian in Fortran (SystemH)
#ifdef _MSC_VER
	void __stdcall HINITHCH4RST(double *coeff, int* diag, int* nmodes);
	void __stdcall HHCH4RST(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void __stdcall INTTOCARTHCH4RST(double* q, double* x);
	void __stdcall KOORDREADHCH4RST();
#endif
}

HCH4_rst::HCH4_rst()
{
}

HCH4_rst::~HCH4_rst()
{
}

void HCH4_rst::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	// Call Hinit here
#ifdef _MSC_VER
	HINITHCH4RST((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	KOORDREADHCH4RST();
#endif
#ifdef __linux__
	hinithch4rst_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	koordreadhch4rst_();
#endif
}

void HCH4_rst::InitOperator()
{
	// set number of parts
	nparts = 164;
	// number of physical modes
	nmodes = 12;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HHCH4RST;
#endif
#ifdef __linux__
	SystemH = hhch4rst_;
#endif
}

Vectord HCH4_rst::IntToCart(const Vectord& q)
{
	Vectord qq(q);
	Vectord Xv(18);
#ifdef _MSC_VER
	INTTOCARTHCH4RST((double*) (&qq(0)), (double*) &Xv(0));
#endif
#ifdef __linux__
	inttocarthch4rst_((double*) (&qq(0)), (double*) &Xv(0));
#endif

	return Xv;
}
