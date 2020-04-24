#include "CH4_radau_J0.h"


extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	// Hamiltonian in Fortran (SystemH)
	void __stdcall HINITCH4RAD(double* coeff, int* diag, int* nmodes);
	void __stdcall HCH4RAD(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void __stdcall INTTOCARTCH4RADAU(double* q, double* x);
#endif

#ifdef __linux__
	void hinitch4rad_(double* coeff, int* diag, int* nmodes);
	void hch4rad_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void inttocartch4radau_(double* q, double* x);
#endif
}


KEO_CH4Rad::KEO_CH4Rad()
{
}


KEO_CH4Rad::~KEO_CH4Rad()
{
}

void KEO_CH4Rad::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	// Call Hinit here
#ifdef _MSC_VER
	HINITCH4RAD((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
#endif
#ifdef __linux__
	hinitch4rad_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
#endif
}

void KEO_CH4Rad::InitOperator()
{
	// set number of parts
	// k zahl=4+34*(N-3)+38*(N-3)*(N-4)+15*(N-3)+1+4
	// for N=6: 384
	// for J=0: k zahl=k zahl-15*(N-3)-4-1
	// which is 50 (for N=6)
	// -> (J=0, N=6) nparts=334, but need 384 for hinit
	// -> (J=0, N=5) nparts=148, but need 183 for hinit (?)
	nparts = 183;
	// number of physical modes
	nmodes = 9;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HCH4RAD;
#endif
#ifdef __linux__
	SystemH = hch4rad_;
#endif
}


Vectord KEO_CH4Rad::IntToCart(const Vectord& q)
{
	Vectord qq(q);
	assert(qq.Dim()==9);
	Vectord Xv(15);
#ifdef _MSC_VER
	INTTOCARTCH4RADAU((double*) (&qq(0)), (double*) &Xv(0));
#endif
#ifdef __linux__
	inttocartch4radau_((double*) (&qq(0)), (double*) &Xv(0));
#endif
	return Xv;
}
