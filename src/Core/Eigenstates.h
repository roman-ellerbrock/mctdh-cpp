//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef EIGENSTATES_H
#define EIGENSTATES_H
#include "Core/HamiltonianRepresentation.h"
#include "Core/IntegratorInterface.h"
#include "TreeClasses/TreeIO.h"

Vectord propagatorEnergies(const Wavefunction& Psi, const Tree& tree, double out) {
	auto S = MatrixTreeFunctions::DotProduct(Psi, Psi, tree);
	const Node& top = tree.TopNode();
	auto propagator_ev = Diagonalize(S[top]).second;
	// E_i=-log(Ev_i)/(2*dt)
	for (int i = 0; i < propagator_ev.Dim(); i++)
		propagator_ev(i) = -log(propagator_ev(i)) / (2. * out);
	return propagator_ev;
}

Vectord Eigenstate(Wavefunction& Psi, const Hamiltonian& H, const Tree& tree) {
	HamiltonianRepresentation hRep(H, tree);
	hRep.build(H, Psi, tree);
	auto Hval = Expectation(hRep, Psi, H, tree);
	auto spec = Diagonalize(Hval);
	auto trafo = spec.first;

	const Node& top = tree.TopNode();
	Psi[top] = multStateArTB(trafo, Psi[top]);
	return spec.second;
}

void Status(const Vectord& eigenvalues, const Vectord& propergatorev,
	ostream& os) {
	pair<string, double> energy("cm", 219474.6313705);

	// Write out Eigenvalues
	int zero = propergatorev.Dim() - 1;
	os << "Eigenvalues (Diagonalizing Hamiltonmatrix"
	   << " / Propagatormatrix)" << endl;
	os << "Ground state :\t" << eigenvalues(0) * energy.second << " " << energy.first;
	os << "\t" << propergatorev(zero) * energy.second << " " << energy.first;
	os << "\t" << eigenvalues(0) << " a.u." << endl;
	for (int i = 1; i < eigenvalues.Dim(); i++) {
		os << i << "\t" << (eigenvalues(i) - eigenvalues(0)) * energy.second << " " << energy.first;
		os << "\t" << (propergatorev(zero - i) - propergatorev(zero)) * energy.second
		   << " " << energy.first;
		os << "\t" << (eigenvalues(i) - eigenvalues(0)) << " a.u." << endl;
	}
}

void Eigenstates(Wavefunction& Psi, const Hamiltonian& H, const Tree& tree) {
	double t = 0.;
	double t_end = 20000.;
	double out = 20.;
	double dt = 1.;
	size_t num_iterations = (t_end - t) / out;

	IntegratorInterface I(H, tree, -QM::im);
	auto energies = Eigenstate(Psi, H, tree);
	Status(energies, energies, cout);
	TreeIO::Output(Psi, tree);

	for (size_t iter = 0; iter < num_iterations; ++iter) {
		// Integrate in imaginary time
		RungeKutta4::Integrate<IntegratorInterface, Wavefunction, double>(t, out, dt, Psi, I);

		auto propagator_energies = propagatorEnergies(Psi, tree, out);

		Orthonormal(Psi, tree);

		// Transform Psi to eigenbasis and get energies
		energies = Eigenstate(Psi, H, tree);

		// I/O
		Status(energies, propagator_energies, cout);
		TreeIO::Output(Psi, tree);
	}
}

#endif //EIGENSTATES_H
