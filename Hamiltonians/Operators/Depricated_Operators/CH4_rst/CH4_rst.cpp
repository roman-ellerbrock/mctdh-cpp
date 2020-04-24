#include "CH4_rst.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	// Hamiltonian in Fortran (SystemH)
	void __stdcall HINITCH4RST(double* coeff, int* diag, int* nmodes);
	// Initialization of masses for the coordinate transformation
	void __stdcall KOORDREADCH4RST();
	// Coordinate Transformation
	void __stdcall INTTOCARTCH4RST(double* q, double* x);
	void __stdcall HCH4RST(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
#ifdef __linux__
	void hinitch4rst_(double* coeff, int* diag, int* nmodes);
	void hch4rst_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void inttocartch4rst_(double* q, double* x);
	void koordreadch4rst_();
#endif
}

CH4_rst::CH4_rst()
{
}

CH4_rst::~CH4_rst()
{
}

void CH4_rst::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
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

void CH4_rst::InitOperator()
{
	// set number of parts
	nparts = 69;
	// number of physical modes
	nmodes = 9;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HCH4RST;
#endif
#ifdef __linux__
	SystemH = hch4rst_;
#endif
}

Vectord CH4_rst::IntToCart(const Vectord& q)
{
	Vectord qq(q);
	Vectord Xv(15);
#ifdef _MSC_VER
	INTTOCARTCH4RST((double*) (&qq(0)), (double*) &Xv(0));
#endif
#ifdef __linux__
	inttocartch4rst_((double*) (&qq(0)), (double*) &Xv(0));
#endif

	return Xv;
}
