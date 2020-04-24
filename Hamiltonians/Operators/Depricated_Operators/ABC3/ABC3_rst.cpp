#include "ABC3_rst.h"

extern "C"
{
	// @TODO: use cmake for this and ensure that it works for all relevant platforms
#ifdef _MSC_VER
	// Hamiltonian in Fortran (SystemH)
	void __stdcall HINITABC3RST(double* coeff, int* diag, int* nmodes);
	// Initialization of masses for the coordinate transformation
	void __stdcall KOORDREADABC3RST(double* mass);
	// Coordinate Transformation
	void __stdcall INTTOCARTABC3RST(double* q, double* x);
	void __stdcall HABC3RST(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
#endif
#ifdef __linux__
	void hinitabc3rst_(double* coeff, int* diag, int* nmodes);
	void habc3rst_(int* mode, int* teil, double* hPsi, double* Psi,
		int* dim, double* matrix, double* trafo, double* ort);
	void inttocartabc3rst_(double* q, double* x);
	void koordreadabc3rst_(double* mass);
#endif
}

ABC3_rst::ABC3_rst(const mctdhBasis &basis)
	:FortranSOP(basis), mass(5)
{
	mass(0)=12.;
	mass(1)=1.0;
	mass(2)=1.0;
	mass(3)=1.0;
	mass(4)=35.446;
}

ABC3_rst::ABC3_rst(const mctdhBasis& basis, Vectord m)
	:FortranSOP(basis), mass(m)
{
	assert(mass.Dim() == 5);
}

ABC3_rst::~ABC3_rst()
{
}

void ABC3_rst::callHinit(Vectorcd & fortrancoeffs, Matrix<int>& diag)
{
	cout << "masses used in operator of ABC_3 molecule:" << endl;
	mass.print();
	// Call Hinit here
#ifdef _MSC_VER
	HINITABC3RST((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	KOORDREADABC3RST(&mass(0));
#endif
#ifdef __linux__
	hinitabc3rst_((double*) (&fortrancoeffs(0)), (int*) &diag(0, 0), &nmodes);
	koordreadabc3rst_(&mass(0));
#endif
}

void ABC3_rst::InitOperator()
{
	// set number of parts
	nparts = 69;
	// number of physical modes
	nmodes = 9;
	// Set Hamiltonian
#ifdef _MSC_VER
	SystemH = HABC3RST;
#endif
#ifdef __linux__
	SystemH = habc3rst_;
#endif
}

/*Vectord ABC3_rst::IntToCart(const Vectord& q)
{
	Vectord qq(q);
	Vectord Xv(15);
#ifdef _MSC_VER
	INTTOCARTABC3RST((double*) (&qq(0)), (double*) &Xv(0));
#endif
#ifdef __linux__
	inttocartabc3rst_((double*) (&qq(0)), (double*) &Xv(0));
#endif

	return Xv;
}
 */
