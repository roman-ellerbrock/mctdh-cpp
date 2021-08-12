//
// Created by Roman Ellerbrock on 6/26/20.
//
#include "Core/tEigenstates.h"
#include "TreeClasses/SpectralDecompositionTree.h"

template<typename T>
Vectord tpropagatorEnergies(const TensorTree<T>& Psi, const Tree& tree, double out) {
	auto S = TreeFunctions::dotProduct(Psi, Psi, tree);
	const Node& top = tree.topNode();
	auto propagator_ev = diagonalize(S[top]).second;
	// E_i=-log(Ev_i)/(2*dt)
	for (int i = 0; i < propagator_ev.dim(); i++)
		propagator_ev(i) = -log(propagator_ev(i)) / (2. * out);
	return propagator_ev;
}

template<typename T>
Vectord tEigenstate(TensorTree<T>& Psi, const SOP<T>& H, const Tree& tree) {

	tHamiltonianRepresentation<T> hRep(H, tree);
	hRep.build(H, Psi, tree);
	auto Hval = Expectation(hRep, Psi, H, tree);
	auto spec = diagonalize(Hval);
	auto trafo = spec.first;

	const Node& top = tree.topNode();
	Psi[top] = multStateArTB(trafo, Psi[top]);
	return spec.second;
}

template<typename T>
void tStatus(const Vectord& eigenvalues, const Vectord& propergatorev,
	const Matrix<T>& S, ostream& os) {
//	pair<string, double> energy("eV", 27.2114);
	pair<string, double> energy("cm", 219474.6313705);

	Vectord s(S.dim1());
	for (size_t i = 0; i < S.dim1(); ++i) {
		for (size_t j = 0; j < S.dim2(); ++j) {
			s(i) += pow(abs(S(i,j)), 2);
		}
	}

	// Write out Eigenvalues
	int zero = propergatorev.dim() - 1;
	std::streamsize prec = os.precision();
	os.precision(12);
	os << "Eigenvalues (Diagonalizing Hamiltonmatrix"
	   << " / Propagatormatrix)" << endl;
	os << "Ground state :\t" << eigenvalues(0) * energy.second << " " << energy.first;
	os << "\t" << propergatorev(zero) * energy.second << " " << energy.first;
	os << "\t" << eigenvalues(0) << " a.u." ;
	os << "\t" << (1.-s(0)) << endl;
	for (int i = 1; i < eigenvalues.dim(); i++) {
		os << i << "\t" << (eigenvalues(i) - eigenvalues(0)) * energy.second << " " << energy.first;
		os << "\t" << (propergatorev(zero - i) - propergatorev(zero)) * energy.second
		   << " " << energy.first;
		os << "\t" << (eigenvalues(i) - eigenvalues(0)) << " a.u.";
		os << "\t(" << (1.-s(i)) << ")" << endl;
	}
	os.precision(prec);
}

template<typename T>
void tEigenstates(tIntegratorVariables<T>& ivar) {
	auto& Psi = *ivar.psi;
	const auto& H = *ivar.sop;
	const Tree& tree = *ivar.tree;
	size_t num_iterations = (ivar.time_end - ivar.time_now) / ivar.out;
	auto eigenvar = ivar;

	cout << "=======================================================" << endl;
	cout << "=============== Eigenstate calculation ================" << endl;
	cout << "=======================================================" << endl;

//	IntegratorInterface I(H, tree, -QM::im);
	tCMFIntegrator<T> cmf(H, tree, -1.);
	auto energies = tEigenstate(Psi, H, tree);
	Matrix<T> s(energies.dim(), energies.dim());
	tStatus(energies, energies, s, cout);
	TreeIO::output2(Psi, tree);

	for (size_t iter = 0; iter < num_iterations; ++iter) {
		TensorTree<T> lastPsi = Psi;

		// Integrate in imaginary time
		eigenvar.time_now=0.;
		eigenvar.time_end =ivar.out -1e-5;
//		RungeKutta4::Integrate<IntegratorInterface, Wavefunction, double>(t, out, dt, Psi, I);
		cmf.Integrate(eigenvar);

		auto propagator_energies = tpropagatorEnergies(Psi, tree, ivar.out);

		orthonormal(Psi, tree);

		// Transform Psi to eigenbasis and get energies
		energies = tEigenstate(Psi, H, tree);
		auto overlap = TreeFunctions::dotProduct(Psi, lastPsi, tree);

		// I/O
		tStatus(energies, propagator_energies, overlap[tree.topNode()], cout);
//		TreeIO::Output(Psi, tree);
	}
}

template void tEigenstates(tIntegratorVariables<double>&);