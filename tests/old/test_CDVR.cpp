//
// Created by Roman Ellerbrock on 3/13/20.
//

#include "DVR/CDVR.h"
#include "Parser/yaml_parser.h"
#include "UnitTest++/UnitTest++.h"
#include "PESs/PESs.h"
#include "Core/HamiltonianRepresentation.h"
#include "Core/Eigenstates.h"

SUITE (CDVR) {
	TEST (coupledHO) {
		string yaml_filename("../examples/test_ExplicitEdgeWavefunction.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		auto Vpot = make_shared<CDVRModelV>(4);
		PotentialOperator V(Vpot, 0, 0);

		CDVR cdvr(tree);
//		CDVR cdvr(Psi, V, tree);

		ExplicitEdgeWavefunction Chi(Psi, tree, true);
	}

	TEST (CDVRHamiltonianRep) {
		string yaml_filename("../examples/cdvr_sl2.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		auto H_ptr = state.hamiltonian_;
		const Hamiltonian H = *H_ptr;

		HamiltonianRepresentation hRep(H, tree);
		hRep.build(H, Psi, tree);

		for (const Node& node : tree) {
			node.info();
			auto Phi = Psi[node];
			auto hPhi = Apply(H, Phi, hRep, node);
			auto hmat = Phi.DotProduct(hPhi);
			hmat *= 219475;
			hmat.print();
		}
	}

	TEST (CDVREigenstate) {
		string yaml_filename("../examples/cdvr_sl2.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		auto H_ptr = state.hamiltonian_;
		const Hamiltonian H = *H_ptr;

//		Eigenstates(Psi, H, tree);
//		cout << "CDVR done.\n" << endl;
	}

	TEST (SOPEigenstate) {

		string yaml_filename("../examples/ho_sl.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		auto H_ptr = state.hamiltonian_;
		const Hamiltonian H = *H_ptr;

//		Eigenstates(Psi, H, tree);
//		cout << "SOP done.\n" << endl;
	}

	TEST (CH3Eigenstate) {

		cout << "ch3 test" << endl;
		string yaml_filename("../examples/ch3.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		auto H_ptr = state.hamiltonian_;
		const Hamiltonian& H = *H_ptr;

/*		const PotentialOperator& Vop = H.V_;
		Vectord Q(6);
		Q(0) = 152.0481;
		Q(1) = 0.95531;
		Q(2) = 0.78539;
		Q(3) = 1.57;
		Q(4) = 1.04719;
		Q(5) = 3.14159;
		double v = Vop.Evaluate(Q, 0);
			CHECK_CLOSE(3.52653e-05, v, 1e-8);
*/
		HamiltonianRepresentation hRep(H, tree);
		cout << "SIZE: " << H.size() << endl;
		hRep.build(H, Psi, tree);
		auto Hval = Expectation(hRep, Psi, H, tree);
		Hval *= 219475.;
			CHECK_CLOSE(6716.77, abs(Hval(0, 0)), 1e-1);

		double t = 0.;
		double t_end = 4000.;
		double out = 100.;
		double dt = 1.0;
		IntegratorVariables ivar(t, t_end, dt, out, 1e-4, 1e-5,
			Psi, H, tree, tree, "out.dat", "in.dat", false);
//		Eigenstates(ivar);
		cout << "CH3 eigenstates done.\n" << endl;
	}
}