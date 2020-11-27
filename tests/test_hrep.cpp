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

	TEST (BuildHmats) {
		string yaml_filename("../examples/coupledHO.ml.yaml");
		auto state = parser::run(yaml_filename);
		HamiltonianRepresentation hrep(*state.hamiltonian_, state.tree_);
		const Hamiltonian& H = *state.hamiltonian_;
		const Wavefunction& Psi = state.wavefunctions_["Psi"];
		const Tree& tree = state.tree_;
		hrep.build(H, Psi, tree);
		for (size_t n = 0; n < 8; ++n) {
			auto& mattree = hrep.hMats_[n];

			auto MPsi = H[n].Apply(Psi, tree);
			auto hmat = TreeFunctions::DotProduct(Psi, MPsi, tree);
			const auto& stree = mattree.Active();
			for (const Node& node : tree) {
				if (!mattree.Active(node)) { continue; }
				TreeFunctions::RepresentLayer(mattree, Psi[node], Psi[node], H[n], node);
				auto residual = Residual(mattree[node], hmat[node]);
					CHECK_CLOSE(0., residual, 1e-7);
			}

			// Check separately
			auto Hexpectation = mattree[tree.TopNode()];
			double value = abs(Hexpectation(0, 0));
			if (n % 2) {
					CHECK_CLOSE(0.5 / 0.018, value, 1e-7);
			} else {
					CHECK_CLOSE(0.25 * 0.018, value, 1e-7);
			}
		}

		for (size_t n = 8; n < 12; ++n) {
			const auto& mattree = hrep.hMats_[n];
			auto MPsi = H[n].Apply(Psi, tree);
			auto hmat = TreeFunctions::DotProduct(Psi, MPsi, tree);
			for (const Node& node : tree) {
				if (!mattree.Active(node)) { continue; }
				auto residual = Residual(mattree[node], hmat[node]);
					CHECK_CLOSE(0., residual, 1e-7);
			}
		}
	}

	TEST (BuildHContractions) {
		string yaml_filename("../examples/coupledho.yaml");
		auto state = parser::run(yaml_filename);
		HamiltonianRepresentation hrep(*state.hamiltonian_, state.tree_);
		const Hamiltonian& H = *state.hamiltonian_;
		const Wavefunction& Psi = state.wavefunctions_["Psi"];
		const Tree& tree = state.tree_;
		hrep.build(H, Psi, tree);

		for (size_t n = 0; n < 8; ++n) {
			const auto& sparse_hcon = hrep.hContractions_[n];
			auto MPsi = H[n].Apply(Psi, tree);
			auto hmat = TreeFunctions::DotProduct(Psi, MPsi, tree);
			auto hcon = TreeFunctions::Contraction(Psi, MPsi, hmat, tree);
			for (const Node& node : tree) {
				if (!sparse_hcon.Active(node)) { continue; }
				auto residual = Residual(sparse_hcon[node], hcon[node]);
					CHECK_CLOSE(0., residual, 1e-7);
			}
		}

		for (size_t n = 8; n < 12; ++n) {
			const auto& sparsehcon = hrep.hContractions_[n];
			auto MPsi = H[n].Apply(Psi, tree);

			auto hmat = TreeFunctions::DotProduct(Psi, MPsi, tree);
			auto hcon = TreeFunctions::Contraction(Psi, MPsi, hmat, tree);
			for (const Node& node : tree) {
				if (!sparsehcon.Active(node)) { continue; }
				auto residual = Residual(sparsehcon[node], hcon[node]);
					CHECK_CLOSE(0., residual, 1e-7);
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
		Hval *= 219475.;
		auto spec = Diagonalize(Hval);
			CHECK_CLOSE(8000., spec.second(0), 5.);
	}
}
