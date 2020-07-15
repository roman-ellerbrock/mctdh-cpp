//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Core/IntegratorInterface.h"
#include "Parser/yaml_parser.h"
#include "Util/RungeKutta4.h"
#include "TreeClasses/TreeIO.h"
#include "Core/Eigenstates.h"
#include "Core/CMFIntegrator.h"

SUITE (IntegratorInterface) {

	/*
	TEST (CoupledHO) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto& Psi = state.wavefunctions_["Psi"];
		const auto& H = *state.hamiltonian_;
		const auto& tree = state.tree_;

		complex<double> phase = -QM::im;
		IntegratorInterface I(H, tree, phase);
		double t = 0.;
		double t_end = 20.;
		double h = 0.2;
		const Node& top = tree.TopNode();

		for (size_t k = 0; k < 100; ++k) {
			RungeKutta4::Integrate<IntegratorInterface, Wavefunction, double>(t, t_end, h, Psi, I);
			auto S = MatrixTreeFunctions::DotProduct(Psi, Psi, tree);
			auto y = Diagonalize(S[top]);
			auto propeigenvalues = y.second;
			// E_i=-log(Ev_i)/(2*dt)
			for (int i = 0; i < propeigenvalues.Dim(); i++)
				propeigenvalues(i) = -log(propeigenvalues(i)) / (2. * t_end);
			propeigenvalues *= QM::cm;

			for (int i = propeigenvalues.Dim() - 1; i >= 0; --i) {
				cout << propeigenvalues(i) << "\t";
			}
			cout << endl;

			Orthogonal(Psi, tree);
			Orthonormal(Psi, tree);

			HamiltonianRepresentation hRep(H, tree);
			hRep.build(H, Psi, tree);
			auto Hval = Expectation(hRep, Psi, H, tree);

			cout << "Energie:\n";
			Hval *= 219475.;
			auto spec = Diagonalize(Hval);
			spec.second.print();
			auto trafo = spec.first;
			Psi[top] = multStateArTB(trafo, Psi[top]);

			hRep.build(H, Psi, tree);
			Hval = Expectation(hRep, Psi, H, tree);
			Hval.print();
			TreeIO::Output(Psi, tree);

			t = 0.;
			if (k == 1) {
				h = 5.;
				t_end = 50.;
			} else if (k == 10) {
				h = 5.;
				t_end = 100.;
			}
		}
	}*/

	TEST(Eigenstates) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto& Psi = state.wavefunctions_["Psi"];
		const auto& H = *state.hamiltonian_;
		const auto& tree = state.tree_;

//		Eigenstates(Psi, H, tree);
	}

	TEST(CMF_Integrator) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		auto& Psi = state.wavefunctions_["Psi"];
		const auto& H = *state.hamiltonian_;
		const auto& tree = state.tree_;

		IntegratorVariables ivar(0., 100., 0.1, 20., 1e-4,
			1e-6, Psi, H, tree, "out.mctdh", "in.mctdh", true);
		complex<double> phase(0., -1.);
		CMFIntegrator cmf(H, tree, phase);
//		cmf.Integrate(ivar, cout);
	}

	TEST(NOCl_Integrator) {

		string yaml_filename("../examples/nocl.yaml");
		auto state = parser::run(yaml_filename);
		auto& Psi = state.wavefunctions_["Psi"];
		const auto& H = *state.hamiltonian_;
		const auto& tree = state.tree_;

		IntegratorVariables ivar(0., 1985.38, 0.1, 41.362, 1e-4,
			1e-4, Psi, H, tree, "out.mctdh", "in.mctdh", true);
		complex<double> phase(1., 0.);
		CMFIntegrator cmf(H, tree, phase);
//		cmf.Integrate(ivar, cout);
	}
}

