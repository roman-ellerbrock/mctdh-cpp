//
// Created by Roman Ellerbrock on 3/13/20.
//

#include "DVR/CDVR.h"
#include "Parser/yaml_parser.h"
#include "UnitTest++/UnitTest++.h"
#include "../Hamiltonians/PESs/PESs.h"
#include "Core/HamiltonianRepresentation.h"
#include "Core/Eigenstates.h"

SUITE(CDVR) {
	TEST(coupledHO) {
		string yaml_filename("../examples/test_ExplicitEdgeWavefunction.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		auto Vpot = make_shared<CDVRModelV>(4);
		PotentialOperator V(Vpot, 0, 0);

		CDVR cdvr(Psi, V, tree);

		ExplicitEdgeWavefunction Chi(Psi, tree, true);

	}

	TEST(CDVRHamiltonianRep) {
		string yaml_filename("../examples/cdvr_sl.yaml");
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

	TEST(CDVREigenstate) {
		string yaml_filename("../examples/cdvr_sl.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		auto H_ptr = state.hamiltonian_;
		const Hamiltonian H = *H_ptr;

		Eigenstates(Psi, H, tree);
		cout << "CDVR done.\n" << endl;
		getchar();
	}

	TEST(SOPEigenstate) {

		string yaml_filename("../examples/ho_sl.yaml");
		auto state = parser::run(yaml_filename);
		auto Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;
		auto H_ptr = state.hamiltonian_;
		const Hamiltonian H = *H_ptr;

		Eigenstates(Psi, H, tree);
		cout << "SOP done.\n" << endl;
		getchar();
	}

}