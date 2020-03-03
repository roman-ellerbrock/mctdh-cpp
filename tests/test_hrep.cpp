//
// Created by Roman Ellerbrock on 3/1/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Core/HamiltonianRepresentation.h"
#include "TreeShape/TreeFactory.h"
#include "Hamiltonians.h"
#include "Parser/yaml_parser.h"

SUITE (HamiltonianRepresentation) {

	TEST (Build) {
		mt19937 gen;
		Tree tree = TreeFactory::BalancedTree(4, 6, 4);
		Hamiltonian H;
		H = CoupledHO(tree);
		HamiltonianRepresentation hrep(H, tree);
		Wavefunction Psi(gen, tree);
		hrep.build(H, Psi, tree);
	}

	TEST (BuildTree) {
		string yaml_filename("../examples/tree.yaml");
		YAML::Node config = YAML::LoadFile(yaml_filename);
		Tree tree = parser::create_tree(config);
	}

	TEST (ReadJobs) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		HamiltonianRepresentation hrep(*state.hamiltonian_, state.tree_);
		hrep.build(*state.hamiltonian_, state.wavefunctions_["Psi"], state.tree_);
		for (size_t n = 0; n < 8; ++n) {
			const auto& mattree = hrep.hMats_[n];
			auto mat = mattree[state.tree_.TopNode()];
			double value = abs(mat(0, 0));
			if (n % 2) {
					CHECK_CLOSE(0.5/0.018, value, 1e-7);
			} else {
					CHECK_CLOSE(0.25*0.018, value, 1e-7);
			}
		}
	}

	TEST (Derivative) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		const auto& H = *state.hamiltonian_;
		const auto& Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		Wavefunction dPsi(tree);
		HamiltonianRepresentation hRep(H, tree);
		double time = 0.;

		Derivative(dPsi, hRep, time, Psi, H, tree);
	}

	TEST (Expectation) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		const auto& H = *state.hamiltonian_;
		const auto& Psi = state.wavefunctions_["Psi"];
		const auto& tree = state.tree_;

		HamiltonianRepresentation hRep(H, tree);
		hRep.build(H, Psi, tree);

		auto Hval = Expectation(hRep, Psi, H, tree);
		cout << "Energie:\n";
		Hval *= 219475.;
		Hval.print();
		auto spec = Diagonalize(Hval);
		spec.second.print();

	}
}
