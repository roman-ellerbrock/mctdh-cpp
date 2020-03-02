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
					CHECK_CLOSE(0.5, value, 1e-7);
			} else {
					CHECK_CLOSE(0.25, value, 1e-7);
			}
		}
	}

	TEST (Derivative) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);

		Wavefunction dPsi(state.tree_);
		HamiltonianRepresentation hRep(*state.hamiltonian_, state.tree_);
		double time = 0.;

		Derivative(dPsi, hRep, time, state.wavefunctions_["Psi"], *state.hamiltonian_, state.tree_);
		dPsi.print(state.tree_);
	}
}
