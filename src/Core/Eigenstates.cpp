//
// Created by Roman Ellerbrock on 6/26/20.
//
#include "Core/Eigenstates.h"
#include "TreeClasses/SpectralDecompositionTree.h"

Vectord propagatorEnergies(const Wavefunction& Psi, const Tree& tree, double out) {
	auto S = TreeFunctions::DotProduct(Psi, Psi, tree);
	const Node& top = tree.TopNode();
	auto propagator_ev = Diagonalize(S[top]).second;
	// E_i=-log(Ev_i)/(2*dt)
	for (int i = 0; i < propagator_ev.Dim(); i++)
		propagator_ev(i) = -log(propagator_ev(i)) / (2. * out);
	return propagator_ev;
}

Vectord Eigenstate(Wavefunction& Psi, const Hamiltonian& H, const Tree& tree,
	const Tree& cdvrtree) {

	HamiltonianRepresentation hRep(H, tree, cdvrtree);
	hRep.build(H, Psi, tree);
	auto Hval = Expectation(hRep, Psi, H, tree);
	auto spec = Diagonalize(Hval);
	auto trafo = spec.first;

	const Node& top = tree.TopNode();
	Psi[top] = multStateArTB(trafo, Psi[top]);
	return spec.second;
}

void Status(const Vectord& eigenvalues, const Vectord& propergatorev,
	const Matrixcd& S, ostream& os) {
//	pair<string, double> energy("eV", 27.2114);
	pair<string, double> energy("cm", 219474.6313705);

	Vectord s(S.Dim1());
	for (size_t i = 0; i < S.Dim1(); ++i) {
		for (size_t j = 0; j < S.Dim2(); ++j) {
			s(i) += pow(abs(S(i,j)), 2);
		}
	}

	// Write out Eigenvalues
	int zero = propergatorev.Dim() - 1;
	os << "Eigenvalues (Diagonalizing Hamiltonmatrix"
	   << " / Propagatormatrix)" << endl;
	os << "Ground state :\t" << eigenvalues(0) * energy.second << " " << energy.first;
	os << "\t" << propergatorev(zero) * energy.second << " " << energy.first;
	os << "\t" << eigenvalues(0) << " a.u." ;
	os << "\t" << (1.-s(0)) << endl;
	for (int i = 1; i < eigenvalues.Dim(); i++) {
		os << i << "\t" << (eigenvalues(i) - eigenvalues(0)) * energy.second << " " << energy.first;
		os << "\t" << (propergatorev(zero - i) - propergatorev(zero)) * energy.second
		   << " " << energy.first;
		os << "\t" << (eigenvalues(i) - eigenvalues(0)) << " a.u.";
		os << "\t(" << (1.-s(i)) << ")" << endl;
	}
}

void Eigenstates(IntegratorVariables& ivar) {
	Wavefunction& Psi = *ivar.psi;
	const Hamiltonian& H = *ivar.h;
	const Tree& tree = *ivar.tree;
	const Tree& cdvrtree = *ivar.cdvrtree;
	size_t num_iterations = (ivar.time_end - ivar.time_now) / ivar.out;
	auto eigenvar = ivar;

	cout << "=======================================================" << endl;
	cout << "=============== Eigenstate calculation ================" << endl;
	cout << "=======================================================" << endl;

//	IntegratorInterface I(H, tree, -QM::im);
	CMFIntegrator cmf(H, tree, cdvrtree, -QM::im);
	auto energies = Eigenstate(Psi, H, tree, cdvrtree);
	Matrixcd s(energies.Dim(), energies.Dim());
	Status(energies, energies, s, cout);
	TreeIO::Output(Psi, tree);

	for (size_t iter = 0; iter < num_iterations; ++iter) {
		Wavefunction lastPsi = Psi;

		// Integrate in imaginary time
		eigenvar.time_now=0.;
		eigenvar.time_end =ivar.out -1e-5;
//		RungeKutta4::Integrate<IntegratorInterface, Wavefunction, double>(t, out, dt, Psi, I);
		cmf.Integrate(eigenvar);

		auto propagator_energies = propagatorEnergies(Psi, tree, ivar.out);

		Orthonormal(Psi, tree);

		// Transform Psi to eigenbasis and get energies
		energies = Eigenstate(Psi, H, tree, cdvrtree);
		auto overlap = TreeFunctions::DotProduct(Psi, lastPsi, tree);

		// I/O
		Status(energies, propagator_energies, overlap[tree.TopNode()], cout);
//		TreeIO::Output(Psi, tree);
	}
}

