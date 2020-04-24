#include "CH4_decoup.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	void __stdcall HINITCH4DECOUP(double* coeff, int* diag, int* nmodes);
#endif
#ifdef __linux__
	void hinitch4decoup_(double* coeff, int* diag, int* nmodes);
#endif
	// Hamiltonian in Fortran (SystemH)
#ifdef _MSC_VER
	void __stdcall HCH4DECOUP(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
#ifdef __linux__
	void hch4decoup_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
	// Coordinate Transformation
#ifdef _MSC_VER
	void __stdcall INTTOCARTCH4DECOUP(double* q, double* x);
#endif
#ifdef __linux__
	void inttocartch4decoup_(double* q, double* x);
#endif
// Initialization of masses for the coordinate transformation
#ifdef _MSC_VER
	void __stdcall KOORDREADCH4DECOUP();
#endif
#ifdef __linux__
	void koordreadch4decoup_();
#endif
}

CH4_decoup::CH4_decoup()
{
}

CH4_decoup::~CH4_decoup()
{
}

void CH4_decoup::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	// Call Hinit here
#ifdef _MSC_VER
	HINITCH4DECOUP((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	KOORDREADCH4DECOUP();
#endif
#ifdef __linux__
	hinitch4decoup_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	koordreadch4decoup_();
#endif
}

void CH4_decoup::InitOperator()
{
	// set number of parts
	nparts = 69;
	// number of physical modes
	nmodes = 5;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HCH4DECOUP;
#endif
#ifdef __linux__
	SystemH = hch4decoup_;
#endif
}

Vectord CH4_decoup::IntToCart(const Vectord& q)
{
	Vectord qq(9);
	qq(0)=q(0);
	qq(1)=0.955317;
	qq(2)=0.785398;
	qq(3)=q(1);
	qq(4)=1.0471976;
	qq(5)=3.14159265;
	qq(6)=q(2);
	qq(7)=q(3);
	qq(8)=q(4);
	Vectord Xv(15);
#ifdef _MSC_VER
	INTTOCARTCH4DECOUP((double*) (&qq(0)), (double*) &Xv(0));
#endif
#ifdef __linux__
	inttocartch4decoup_((double*) (&qq(0)), (double*) &Xv(0));
#endif

	return Xv;
}
